#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "GeometryShop.H"
#include "BaseIF.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "CH_Timer.H"
#include "memusage.H"
#include "memtrack.H"
#include "EBEllipticLoadBalance.H"
#include "IntVectSet.H"

using std::cerr;
using std::cout;
using std::endl;

/**
 * Direct Array Interface for SDF data
 * Implements BaseIF to directly use numpy array data
 */
class DirectArrayIF: public BaseIF
{
public:
  DirectArrayIF(const std::vector<Real>& a_data,
                const IntVect& a_resolution,
                const RealVect& a_dx,
                const RealVect& a_origin,
                const Real& a_levelValue,
                const bool& a_inside);

  virtual ~DirectArrayIF();

  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  std::vector<Real> m_data;
  IntVect m_resolution;
  RealVect m_dx;
  RealVect m_origin;
  Real m_levelValue;
  bool m_inside;

  // Helper function for trilinear interpolation
  Real interpolateValue(const RealVect& a_point) const;
  int getIndex(int i, int j, int k) const;
};

DirectArrayIF::DirectArrayIF(const std::vector<Real>& a_data,
                             const IntVect& a_resolution,
                             const RealVect& a_dx,
                             const RealVect& a_origin,
                             const Real& a_levelValue,
                             const bool& a_inside)
{
  m_data = a_data;
  m_resolution = a_resolution;
  m_dx = a_dx;
  m_origin = a_origin;
  m_levelValue = a_levelValue;
  m_inside = a_inside;
}

DirectArrayIF::~DirectArrayIF()
{
}

int DirectArrayIF::getIndex(int i, int j, int k) const
{
  // Clamp indices to bounds
  i = std::max(0, std::min(i, m_resolution[0] - 1));
  j = std::max(0, std::min(j, m_resolution[1] - 1));
  k = std::max(0, std::min(k, m_resolution[2] - 1));
  
  return i + j * m_resolution[0] + k * m_resolution[0] * m_resolution[1];
}

Real DirectArrayIF::interpolateValue(const RealVect& a_point) const
{
  // Convert physical coordinates to grid coordinates
  RealVect gridCoord = (a_point - m_origin) / m_dx;
  
  // Get integer grid indices
  int i0 = (int)floor(gridCoord[0]);
  int j0 = (int)floor(gridCoord[1]);
  int k0 = (int)floor(gridCoord[2]);
  
  int i1 = i0 + 1;
  int j1 = j0 + 1;
  int k1 = k0 + 1;
  
  // Get fractional parts
  Real fx = gridCoord[0] - i0;
  Real fy = gridCoord[1] - j0;
  Real fz = gridCoord[2] - k0;
  
  // Trilinear interpolation
  Real v000 = m_data[getIndex(i0, j0, k0)];
  Real v001 = m_data[getIndex(i0, j0, k1)];
  Real v010 = m_data[getIndex(i0, j1, k0)];
  Real v011 = m_data[getIndex(i0, j1, k1)];
  Real v100 = m_data[getIndex(i1, j0, k0)];
  Real v101 = m_data[getIndex(i1, j0, k1)];
  Real v110 = m_data[getIndex(i1, j1, k0)];
  Real v111 = m_data[getIndex(i1, j1, k1)];
  
  Real v00 = v000 * (1 - fx) + v100 * fx;
  Real v01 = v001 * (1 - fx) + v101 * fx;
  Real v10 = v010 * (1 - fx) + v110 * fx;
  Real v11 = v011 * (1 - fx) + v111 * fx;
  
  Real v0 = v00 * (1 - fy) + v10 * fy;
  Real v1 = v01 * (1 - fy) + v11 * fy;
  
  Real value = v0 * (1 - fz) + v1 * fz;
  
  return value;
}

Real DirectArrayIF::value(const RealVect& a_point) const
{
  Real sdfValue = interpolateValue(a_point);
  Real distance = sdfValue - m_levelValue;
  
  if (m_inside)
  {
    return -distance;
  }
  else
  {
    return distance;
  }
}

BaseIF* DirectArrayIF::newImplicitFunction() const
{
  DirectArrayIF* newIF = new DirectArrayIF(m_data, m_resolution, m_dx, m_origin,
                                           m_levelValue, m_inside);
  return static_cast<BaseIF*>(newIF);
}

/**
 * NPZ/Binary File Reader
 *
 * Supports two formats:
 * 1. Python .npz (ZIP archive containing .npy files)
 * 2. Simple binary format (header + data)
 */
class NPZReader
{
public:
  static std::vector<Real> readNPZ(const std::string& filename, IntVect& resolution);
  static std::vector<Real> readBinaryFormat(const std::string& filename, IntVect& resolution);
  static std::vector<Real> readNumpyNPZ(const std::string& filename, IntVect& resolution);

private:
  static bool isZipFile(const std::string& filename);
  static std::vector<Real> parseNPYArray(const char* npy_data, size_t npy_size, IntVect& resolution);
};

bool NPZReader::isZipFile(const std::string& filename)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) return false;

  char magic[4];
  file.read(magic, 4);
  file.close();

  // ZIP files start with "PK\x03\x04"
  return (magic[0] == 'P' && magic[1] == 'K' && magic[2] == 0x03 && magic[3] == 0x04);
}

std::vector<Real> NPZReader::parseNPYArray(const char* npy_data, size_t npy_size, IntVect& resolution)
{
  // NPY format: magic (6 bytes) + version (2 bytes) + header_len + header_dict + data
  if (npy_size < 10 || memcmp(npy_data, "\x93NUMPY", 6) != 0)
  {
    MayDay::Error("Invalid NPY format");
  }

  unsigned char major = npy_data[6];
  unsigned char minor = npy_data[7];

  int header_len;
  if (major == 1)
  {
    header_len = *(unsigned short*)(npy_data + 8);
  }
  else if (major == 2 || major == 3)
  {
    header_len = *(unsigned int*)(npy_data + 8);
  }
  else
  {
    MayDay::Error("Unsupported NPY version");
  }

  int header_offset = (major == 1) ? 10 : 12;
  std::string header(npy_data + header_offset, header_len);

  pout() << "NPY header: " << header << endl;

  // Parse shape from header (looking for 'shape': (dim0, dim1, dim2))
  size_t shape_pos = header.find("'shape'");
  if (shape_pos == std::string::npos)
  {
    MayDay::Error("Could not find shape in NPY header");
  }

  size_t paren_start = header.find('(', shape_pos);
  size_t paren_end = header.find(')', paren_start);
  std::string shape_str = header.substr(paren_start + 1, paren_end - paren_start - 1);

  // Parse dimensions
  int dims[3];
  int dim_count = 0;
  size_t pos = 0;
  while (pos < shape_str.length() && dim_count < 3)
  {
    while (pos < shape_str.length() && !isdigit(shape_str[pos])) pos++;
    if (pos >= shape_str.length()) break;

    dims[dim_count++] = atoi(shape_str.c_str() + pos);
    while (pos < shape_str.length() && isdigit(shape_str[pos])) pos++;
  }

  // Handle both 1D flattened and 3D arrays
  int data_size;
  if (dim_count == 1)
  {
    // Flattened array - compute cube root to get resolution
    data_size = dims[0];
    double cube_root = pow(data_size, 1.0/3.0);
    int res = (int)(cube_root + 0.5); // Round to nearest int

    if (res * res * res != data_size)
    {
      pout() << "Data size " << data_size << " is not a perfect cube" << endl;
      MayDay::Error("Flattened array must have size that is a perfect cube");
    }

    pout() << "Reshaping 1D array (" << data_size << ",) to 3D (" << res << ", " << res << ", " << res << ")" << endl;
    dims[0] = dims[1] = dims[2] = res;
  }
  else if (dim_count == 3)
  {
    if (dims[0] != dims[1] || dims[1] != dims[2])
    {
      MayDay::Error("Non-cubic data not supported");
    }
    data_size = dims[0] * dims[1] * dims[2];
  }
  else
  {
    MayDay::Error("Expected 1D (flattened) or 3D array");
  }

  resolution = IntVect(D_DECL(dims[0], dims[1], dims[2]));

  pout() << "Array shape: (" << dims[0] << ", " << dims[1] << ", " << dims[2] << ")" << endl;

  // Determine data type from header
  bool is_float32 = (header.find("'<f4'") != std::string::npos || header.find("'<f'") != std::string::npos);
  bool is_float64 = (header.find("'<f8'") != std::string::npos);

  if (!is_float32 && !is_float64)
  {
    // Try to auto-detect from data size
    const char* data_ptr = npy_data + header_offset + header_len;
    size_t data_bytes = npy_size - (header_offset + header_len);

    if (data_bytes >= data_size * sizeof(double))
    {
      is_float64 = true;
      pout() << "Auto-detected float64 from data size" << endl;
    }
    else if (data_bytes >= data_size * sizeof(float))
    {
      is_float32 = true;
      pout() << "Auto-detected float32 from data size" << endl;
    }
    else
    {
      MayDay::Error("Cannot determine data type or insufficient data");
    }
  }

  pout() << "Data type: " << (is_float32 ? "float32" : "float64") << endl;

  const char* data_ptr = npy_data + header_offset + header_len;
  std::vector<Real> data(data_size);

  if (is_float32)
  {
    const float* float_data = reinterpret_cast<const float*>(data_ptr);
    // For 1D flattened arrays, data is already in row-major (C) order
    // For 3D arrays from Python, transpose from (z,y,x) to (x,y,z)
    if (dim_count == 1)
    {
      // Already in correct C order
      for (int i = 0; i < data_size; i++)
      {
        data[i] = static_cast<Real>(float_data[i]);
      }
    }
    else
    {
      // Transpose 3D array
      for (int k = 0; k < dims[2]; k++)
      {
        for (int j = 0; j < dims[1]; j++)
        {
          for (int i = 0; i < dims[0]; i++)
          {
            int py_idx = k * dims[0] * dims[1] + j * dims[0] + i;
            int cpp_idx = i + j * dims[0] + k * dims[0] * dims[1];
            data[cpp_idx] = static_cast<Real>(float_data[py_idx]);
          }
        }
      }
    }
  }
  else // float64
  {
    const double* double_data = reinterpret_cast<const double*>(data_ptr);
    if (dim_count == 1)
    {
      // Already in correct C order
      for (int i = 0; i < data_size; i++)
      {
        data[i] = static_cast<Real>(double_data[i]);
      }
    }
    else
    {
      // Transpose 3D array
      for (int k = 0; k < dims[2]; k++)
      {
        for (int j = 0; j < dims[1]; j++)
        {
          for (int i = 0; i < dims[0]; i++)
          {
            int py_idx = k * dims[0] * dims[1] + j * dims[0] + i;
            int cpp_idx = i + j * dims[0] + k * dims[0] * dims[1];
            data[cpp_idx] = static_cast<Real>(double_data[py_idx]);
          }
        }
      }
    }
  }

  return data;
}

std::vector<Real> NPZReader::readNumpyNPZ(const std::string& filename, IntVect& resolution)
{
  pout() << "Reading NumPy .npz file (ZIP archive format)..." << endl;

  // Create unique temporary file name for each MPI rank
  int rank = 0;
#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
#endif
  std::string temp_npy = "/tmp/temp_data_rank" + std::to_string(rank) + ".npy";

  // Use system command to extract the .npy file from the .npz archive
  // Try "sdf.npy" first, then "data.npy" if that fails
  std::string cmd = "unzip -p \"" + filename + "\" sdf.npy > " + temp_npy + " 2>/dev/null";

  int ret = system(cmd.c_str());
  if (ret != 0)
  {
    // Try "data.npy" instead
    pout() << "  'sdf.npy' not found, trying 'data.npy'..." << endl;
    cmd = "unzip -p \"" + filename + "\" data.npy > " + temp_npy + " 2>/dev/null";
    ret = system(cmd.c_str());
    if (ret != 0)
    {
      MayDay::Error("Failed to extract sdf.npy or data.npy from .npz archive. Is 'unzip' available?");
    }
  }

  // Read the extracted .npy file
  std::ifstream file(temp_npy, std::ios::binary);
  if (!file.is_open())
  {
    MayDay::Error("Cannot open extracted .npy file");
  }

  // Get file size
  file.seekg(0, std::ios::end);
  size_t file_size = file.tellg();
  file.seekg(0, std::ios::beg);

  // Read entire file into memory
  std::vector<char> npy_buffer(file_size);
  file.read(npy_buffer.data(), file_size);
  file.close();

  // Clean up temp file
  remove(temp_npy.c_str());

  // Parse the NPY data
  std::vector<Real> data = parseNPYArray(npy_buffer.data(), file_size, resolution);

  // Print statistics
  Real min_val = *std::min_element(data.begin(), data.end());
  Real max_val = *std::max_element(data.begin(), data.end());
  pout() << "Data range: [" << min_val << ", " << max_val << "]" << endl;

  return data;
}

std::vector<Real> NPZReader::readBinaryFormat(const std::string& filename, IntVect& resolution)
{
  pout() << "Reading simple binary format..." << endl;

  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open())
  {
    MayDay::Error("Cannot open binary data file");
  }

  // Read header: resolution (3 ints) + data_size (1 int)
  int res_x, res_y, res_z, data_size;
  file.read(reinterpret_cast<char*>(&res_x), sizeof(int));
  file.read(reinterpret_cast<char*>(&res_y), sizeof(int));
  file.read(reinterpret_cast<char*>(&res_z), sizeof(int));
  file.read(reinterpret_cast<char*>(&data_size), sizeof(int));

  if (res_x != res_y || res_y != res_z)
  {
    MayDay::Error("Non-cubic data not supported");
  }

  resolution = IntVect(D_DECL(res_x, res_y, res_z));

  pout() << "Reading binary data: " << res_x << "^3 = " << data_size << " points" << endl;

  // Read data as double precision
  std::vector<double> temp_data(data_size);
  file.read(reinterpret_cast<char*>(temp_data.data()), data_size * sizeof(double));

  if (!file.good())
  {
    MayDay::Error("Error reading binary data file");
  }

  file.close();

  // Convert to Real (which might be float or double)
  std::vector<Real> data(data_size);
  for (int i = 0; i < data_size; i++)
  {
    data[i] = static_cast<Real>(temp_data[i]);
  }

  // Print statistics
  Real min_val = *std::min_element(data.begin(), data.end());
  Real max_val = *std::max_element(data.begin(), data.end());
  pout() << "Data range: [" << min_val << ", " << max_val << "]" << endl;

  return data;
}

std::vector<Real> NPZReader::readNPZ(const std::string& filename, IntVect& resolution)
{
  // Auto-detect file format
  if (isZipFile(filename))
  {
    return readNumpyNPZ(filename, resolution);
  }
  else
  {
    return readBinaryFormat(filename, resolution);
  }
}

/**
 * EBIS Creation with Adaptive Mesh Refinement (HDF5ToEBIS-style)
 * Uses geometry-aware tagging and BRMeshRefine
 */
void createEBISFromArrayAdaptive(const std::vector<Real>& data,
                                 const IntVect& resolution,
                                 const std::string& outputFile,
                                 Real domainLength = 2.0,
                                 Real levelValue = 0.0,
                                 bool insideRegular = false,
                                 int maxGridSize = 64,
                                 int maxLevel = 0,
                                 int refRatio = 2,
                                 int baseNCells = 16,
                                 Real fillRatio = 0.7,
                                 int blockFactor = 4,
                                 int bufferSize = 1,
                                 int ebTagBuffer = 0)
{
  CH_TIME("createEBISFromArrayAdaptive");

  // Setup domain parameters
  RealVect origin = RealVect::Zero;
  RealVect sdfDx = RealVect::Unit * (domainLength / resolution[0]);

  // Set up base domain
  IntVect lo = IntVect::Zero;
  IntVect hi = IntVect(D_DECL(baseNCells-1, baseNCells-1, baseNCells-1));
  Box baseDomainBox(lo, hi);
  ProblemDomain baseDomain(baseDomainBox);
  RealVect baseDx = RealVect::Unit * (domainLength / baseNCells);

  // Calculate finest domain
  int refinementFactor = 1;
  for (int lev = 0; lev < maxLevel; lev++)
  {
    refinementFactor *= refRatio;
  }
  int finestNCells = baseNCells * refinementFactor;
  IntVect finestHi = IntVect(D_DECL(finestNCells-1, finestNCells-1, finestNCells-1));
  Box finestBox(IntVect::Zero, finestHi);
  ProblemDomain finestDomain(finestBox);
  RealVect finestDx = RealVect::Unit * (domainLength / finestNCells);

  pout() << "Creating EBIS with Adaptive AMR (HDF5ToEBIS-style):" << endl;
  pout() << "  SDF Resolution: " << resolution << endl;
  pout() << "  SDF dx: " << sdfDx << endl;
  pout() << "  Base resolution: " << baseNCells << endl;
  pout() << "  Max Level: " << maxLevel << endl;
  pout() << "  Refinement factor: " << refinementFactor << endl;
  pout() << "  Finest resolution: " << finestNCells << endl;
  pout() << "  Finest dx: " << finestDx << endl;
  pout() << "  Fill ratio: " << fillRatio << endl;
  pout() << "  Block factor: " << blockFactor << endl;
  pout() << "  Buffer size: " << bufferSize << endl;
  pout() << "  EB tag buffer: " << ebTagBuffer << endl;

  // Create the implicit function from array data
  DirectArrayIF arrayIF(data, resolution, sdfDx, origin, levelValue, insideRegular);

  // Create geometry shop
  int verbosity = 1;
  GeometryShop geometryShop(arrayIF, verbosity, finestDx);

  // First, define EBIS at finest level to get geometry
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int ebMaxCoarsen = -1;
  ebisPtr->define(finestDomain, origin, finestDx[0], geometryShop, maxGridSize, ebMaxCoarsen);

  pout() << "EBIS defined at finest level, now building adaptive grids..." << endl;

  // Now build the adaptive grids using BRMeshRefine
  // Split up coarsest domain by max box size
  Vector<Box> boxesCoarsest;
  domainSplit(baseDomain, boxesCoarsest, maxGridSize, blockFactor);

  Vector<int> procsCoarsest;
  LoadBalance(procsCoarsest, boxesCoarsest);
  DisjointBoxLayout dblCoarsest(boxesCoarsest, procsCoarsest);

  // Make EBISLayout at coarsest level
  EBISLayout ebislCoarsest;
  ebisPtr->fillEBISLayout(ebislCoarsest, dblCoarsest, baseDomain, 0);

  // Tag all irregular cells at coarsest level
  pout() << "Tagging irregular cells at coarsest level..." << endl;
  IntVectSet tagsCoarsestLocal;
  int totalIrregCells = 0;

  for (DataIterator dit = dblCoarsest.dataIterator(); dit.ok(); ++dit)
  {
    const EBISBox& ebisBox = ebislCoarsest[dit()];
    const Box& box = dblCoarsest.get(dit());
    IntVectSet irregIVS = ebisBox.getIrregIVS(box);
    tagsCoarsestLocal |= irregIVS;
    totalIrregCells += irregIVS.numPts();
  }

  pout() << "Found " << totalIrregCells << " irregular cells (cells cut by geometry interface)" << endl;

  // Show spatial extent of tagged cells
  if (tagsCoarsestLocal.numPts() > 0)
  {
    Box tagBox = tagsCoarsestLocal.minBox();
    pout() << "Tagged cell bounding box: " << tagBox << endl;
    pout() << "  Physical extent: ";
    pout() << "X=[" << tagBox.smallEnd(0)*baseDx[0] << "," << (tagBox.bigEnd(0)+1)*baseDx[0] << "], ";
    pout() << "Y=[" << tagBox.smallEnd(1)*baseDx[1] << "," << (tagBox.bigEnd(1)+1)*baseDx[1] << "], ";
    pout() << "Z=[" << tagBox.smallEnd(2)*baseDx[2] << "," << (tagBox.bigEnd(2)+1)*baseDx[2] << "]" << endl;
  }

  // Buffer the irregular cells
  if (ebTagBuffer > 0)
  {
    pout() << "Buffering tags by " << ebTagBuffer << " cells" << endl;
    tagsCoarsestLocal.grow(ebTagBuffer);
  }

  pout() << "Total tagged cells (after buffering): " << tagsCoarsestLocal.numPts() << endl;

  // Generate vector of refinement ratios
  Vector<int> refRatios(maxLevel + 1, refRatio);

  // Use BRMeshRefine to generate adaptive grids
  BRMeshRefine gridder(baseDomain, refRatios, fillRatio, blockFactor, bufferSize, maxGridSize);

  Vector<Vector<Box> > newMeshes(maxLevel + 1);
  Vector<Vector<Box> > oldMeshes(maxLevel + 1);
  oldMeshes[0] = boxesCoarsest;

  for (int ilev = 1; ilev <= maxLevel; ilev++)
  {
    oldMeshes[ilev] = Vector<Box>(1, refine(oldMeshes[ilev-1][0], refRatio));
  }

  int baseLevel = 0;
  pout() << "Running BRMeshRefine to create adaptive grids..." << endl;
  gridder.regrid(newMeshes, tagsCoarsestLocal, baseLevel, maxLevel, oldMeshes);

  // Print mesh statistics
  long long totalPoints = 0;
  long long totalBoxes = 0;
  for (int ilev = 0; ilev <= maxLevel; ilev++)
  {
    long long pointsThisLevel = 0;
    for (int ibox = 0; ibox < newMeshes[ilev].size(); ibox++)
    {
      pointsThisLevel += newMeshes[ilev][ibox].numPts();
    }
    totalPoints += pointsThisLevel;
    totalBoxes += newMeshes[ilev].size();
    pout() << "Level " << ilev << ": " << newMeshes[ilev].size() << " boxes, "
           << pointsThisLevel << " points" << endl;
  }
  pout() << "Total: " << totalBoxes << " boxes, " << totalPoints << " points" << endl;

  // Write EBIS to file
#ifdef CH_USE_HDF5
  HDF5Handle geomHandle(outputFile, HDF5Handle::CREATE);
  pout() << "Writing EBIS to file: " << outputFile << endl;
  ebisPtr->write(geomHandle);
  geomHandle.close();
#else
  pout() << "HDF5 not enabled - cannot write EBIS file" << endl;
#endif

  pout() << "Adaptive EBIS creation completed!" << endl;
}

/**
 * Main EBIS Creation Function with AMR support (original uniform approach)
 */
void createEBISFromArray(const std::vector<Real>& data,
                        const IntVect& resolution,
                        const std::string& outputFile,
                        Real domainLength = 2.0,
                        Real levelValue = 0.0,
                        bool insideRegular = false,
                        int maxGridSize = 64,
                        int maxLevel = 0,
                        int refRatio = 2,
                        int baseNCells = 16)
{
  CH_TIME("createEBISFromArray");

  // The NPZ data resolution represents the implicit function sampling,
  // which we'll use for the GeometryShop
  RealVect origin = RealVect::Zero;
  RealVect sdfDx = RealVect::Unit * (domainLength / resolution[0]);

  // Calculate the finest level resolution
  // The three-stage method defines EBIS at max_level resolution
  int refinementFactor = 1;
  for (int lev = 0; lev < maxLevel; lev++)
  {
    refinementFactor *= refRatio;
  }

  int finestNCells = baseNCells * refinementFactor;

  // Set up the computational domain at the finest level
  IntVect lo = IntVect::Zero;
  IntVect hi = IntVect(D_DECL(finestNCells-1, finestNCells-1, finestNCells-1));
  Box finestBox(lo, hi);
  ProblemDomain finestDomain(finestBox);

  // Finest dx based on the computational grid, not the SDF grid
  RealVect finestDx = RealVect::Unit * (domainLength / finestNCells);

  pout() << "Creating EBIS with AMR:" << endl;
  pout() << "  SDF Resolution: " << resolution << endl;
  pout() << "  SDF dx: " << sdfDx << endl;
  pout() << "  Base resolution: " << baseNCells << endl;
  pout() << "  Max Level: " << maxLevel << endl;
  pout() << "  Refinement factor: " << refinementFactor << endl;
  pout() << "  Finest resolution: " << finestNCells << endl;
  pout() << "  Finest Domain: " << finestDomain << endl;
  pout() << "  Finest dx: " << finestDx << endl;
  pout() << "  Level value: " << levelValue << endl;
  pout() << "  Inside regular: " << insideRegular << endl;

  // Create the implicit function from array data
  DirectArrayIF arrayIF(data, resolution, sdfDx, origin, levelValue, insideRegular);

  // Create geometry shop - use the SDF resolution for geometry queries
  int verbosity = 1;
  GeometryShop geometryShop(arrayIF, verbosity, sdfDx);

  // Create EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int ebMaxCoarsen = -1;

  // Define at the finest level to match the three-stage method
  ebisPtr->define(finestDomain, origin, finestDx[0], geometryShop, maxGridSize, ebMaxCoarsen);

  // Write EBIS to file
#ifdef CH_USE_HDF5
  HDF5Handle geomHandle(outputFile, HDF5Handle::CREATE);
  pout() << "Writing EBIS to file: " << outputFile << endl;
  ebisPtr->write(geomHandle);
  geomHandle.close();
#else
  pout() << "HDF5 not enabled - cannot write EBIS file" << endl;
#endif

  pout() << "EBIS creation completed!" << endl;
}

/**
 * Main function
 */
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  {
    CH_TIMERS("npz_to_ebis_direct");
    CH_TIMER("total", t1);
    
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    CH_START(t1);
    
    if (argc < 3)
    {
      cerr << "Usage: " << argv[0] << " <input.npz> <output.ebis.hdf5> [level_value] [inside_regular] [max_level] [use_adaptive]" << endl;
      cerr << "Example: " << argv[0] << " cube.npz cube.ebis.hdf5 0.0 0 3 1" << endl;
      cerr << "  max_level: number of AMR levels (0 = single level, 3 = 4 levels total)" << endl;
      cerr << "  use_adaptive: 0 = uniform refinement (default), 1 = adaptive (HDF5ToEBIS-style)" << endl;
      exit(1);
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    Real levelValue = (argc > 3) ? atof(argv[3]) : 0.0;
    bool insideRegular = (argc > 4) ? (atoi(argv[4]) != 0) : false;
    int maxLevel = (argc > 5) ? atoi(argv[5]) : 0;
    bool useAdaptive = (argc > 6) ? (atoi(argv[6]) != 0) : false;

    try
    {
      // Read NPZ data
      IntVect resolution;
      std::vector<Real> data = NPZReader::readNPZ(inputFile, resolution);

      pout() << "Loaded " << data.size() << " data points from " << inputFile << endl;

      if (useAdaptive)
      {
        pout() << "Using ADAPTIVE mesh refinement (HDF5ToEBIS-style)" << endl;
        // Use adaptive mesh refinement with tagging
        // Use baseNCells=32 for base resolution with user-provided max_level
        createEBISFromArrayAdaptive(data, resolution, outputFile, 2.0, levelValue, insideRegular,
                                   32, maxLevel, 2, 32, 0.7, 4, 1, 0);
      }
      else
      {
        pout() << "Using UNIFORM mesh refinement" << endl;
        // Create EBIS directly from array data with uniform AMR
        createEBISFromArray(data, resolution, outputFile, 2.0, levelValue, insideRegular, 32, maxLevel, 2, 32);
      }

    }
    catch (const std::exception& e)
    {
      cerr << "Error: " << e.what() << endl;
      return 1;
    }
    
    CH_STOP(t1);
  }

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif

  return 0;
} 