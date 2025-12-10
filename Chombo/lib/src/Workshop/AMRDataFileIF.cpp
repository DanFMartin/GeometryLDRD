#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRDataFileIF.H"
#include "AMRIO.H"
#include "QuadCFInterp.H"

#include "NamespaceHeader.H"

using std::cin;

AMRDataFileIF::AMRDataFileIF(const char* const           a_filename,
                             const Real&                 a_value,
                             const bool&                 a_inside,
                             const bool&                 a_useCubicInterp)
{

  // Read all the data from the file
  ReadData(m_noDataValue,a_filename);

  // Save some state
  m_value = a_value;
  m_inside = a_inside;
  m_useCubicInterp = a_useCubicInterp;

  // Make the unit cube IntVectSet
  MakeCorners();
}



AMRDataFileIF::AMRDataFileIF(const AMRDataFileIF& a_inputIF)
{
  // Save a refcounted pointer to the data
  m_real_data   = a_inputIF.m_real_data;
  m_char_data   = a_inputIF.m_char_data;
  m_noDataValue = a_inputIF.m_noDataValue;

  // Copy all the other data
  m_baseDomain = a_inputIF.m_baseDomain;

  m_vectDx = a_inputIF.m_vectDx;
  m_origin = a_inputIF.m_origin;
  m_refRatio = a_inputIF.m_refRatio;
  
  m_value = a_inputIF.m_value;
  m_inside = a_inputIF.m_inside;
  m_useCubicInterp = a_inputIF.m_useCubicInterp;

  // Make the unit cube IntVectSet
  MakeCorners();
}

AMRDataFileIF::AMRDataFileIF(const RefCountedPtr<Vector<LevelData<FArrayBox>* > > a_real_data,
                             const RefCountedPtr<BaseFab<unsigned char> > a_char_data,
                             const Real&                                  a_noDataValue,
                             const Box&                                   a_baseDomain,
                             const Vector<RealVect>&                      a_vectDx,
                             const RealVect&                              a_origin,
                             const Vector<int>&                           a_refRatio,
                             const Real&                                  a_value,
                             const bool&                                  a_inside,
                             const bool&                                  a_useCubicInterp)
{
  // Save a refcounted pointer to the data
  m_real_data   = a_real_data;
  m_char_data   = a_char_data;
  m_noDataValue = a_noDataValue;
  m_numLevels = m_real_data->size();
  
  // Copy all the other data
  m_baseDomain = a_baseDomain;

  m_vectDx = a_vectDx;
  m_origin = a_origin;
  m_refRatio = a_refRatio;
  
  m_value = a_value;
  m_inside = a_inside;
  m_useCubicInterp = a_useCubicInterp;

  // Make the unit cube IntVectSet
  MakeCorners();
}

AMRDataFileIF::~AMRDataFileIF()
{
}

void AMRDataFileIF::GetHeader(Box&  a_baseDomain,
                              Vector<RealVect>& a_vectDx,
                              RealVect& a_origin) const
{
  // Copy header information over
  a_baseDomain = m_baseDomain;
  a_vectDx = m_vectDx;
  a_origin = m_origin;
}

void AMRDataFileIF::GetParams(Real& a_value,
                              bool& a_inside,
                              bool& a_useCubicInterp) const
{
  // Copy parameter information over
  a_value = m_value;
  a_inside = m_inside;
  a_useCubicInterp = m_useCubicInterp;
}

void AMRDataFileIF::SetParams(const Real& a_value,
                              const bool& a_inside,
                              const bool& a_useCubicInterp)
{
  // Set parameter information
  m_value = a_value;
  m_inside = a_inside;
  m_useCubicInterp = a_useCubicInterp;
}

void AMRDataFileIF::SetNoDataValue(const Real& a_value)
{
  // value to use when we are outside
  m_noDataValue = a_value;
}

Real AMRDataFileIF::value(const RealVect& a_point) const
{
  IndexTM<Real,GLOBALDIM> point;

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    point[idir] = a_point[idir];
  }

  for (int idir = SpaceDim; idir < GLOBALDIM; idir++)
  {
    point[idir] = 0.0;
  }

  return value(point);
}

Real AMRDataFileIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  return value(IndexTM<int,GLOBALDIM>::Zero,a_point);
}

Real AMRDataFileIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                          const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Real retval;

  // drop down through the AMR levels and locate valid region which contains this point

  // rachet domain box up to finest level
  Box levelDomainBox = m_baseDomain;
  for (int level=1; level<m_numLevels; level++)
    {
      levelDomainBox.refine(m_refRatio[level-1]);
    }

  
  bool found = false;
  int lev = m_numLevels -1;
  
  while (!found)
    {
      
      // The index in the data corresponding to a_point
      IntVect loCorner;
      // The fraction of a cell from the low index to a_point
      RealVect fraction;
      
      // Determine the low index (which cell the point is in) and the fraction
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real index;
          
          // The floating point "index" of a_point in the data
          index = (a_point[idir] - m_origin[idir]) / m_vectDx[lev][idir];
          
          // The integer portion of that "index"
          loCorner[idir] = (int)floor(index);
        
          // Make sure the low corner in one inside the high end of the data
          // by flagging if we're on the high edge of the base domain on level 0
          // might need to finess this a bit if 
          if ((loCorner[idir] == (levelDomainBox.bigEnd(idir)+1)))
            {
              loCorner[idir] -= 1;
            }

          
          // The fraction portion of the "index"
          fraction[idir] = index - loCorner[idir];
        }
      
      // DFM (3/17/25) I don't understand why this is here
#if 0
      for (int idir = SpaceDim; idir < GLOBALDIM; idir++)
        {
          Real index;
          
          // The floating point "index" of a_point in the data
          index = (a_point[idir] - m_origin[idir]) / m_spacing[idir];
          
          // The integer portion of that "index"
          int loCorner = (int)floor(index);
          
          // The fraction portion of the "index"
          fraction[idir] = index - loCorner;
        }
#endif
    
      // now see if this cell is in the valid domain for this level
      const LevelData<FArrayBox>& levelLDF = *(m_real_data->operator[](lev));
      const DisjointBoxLayout& levelGrids = levelLDF.getBoxes();
      DataIterator dit =levelGrids.dataIterator();
      for (dit.begin(); (dit.ok() && !found); ++dit)
        {
          const Box& validBox = levelGrids[dit];
          const Box dataBox = levelLDF[dit].box();
          
          // does this valid box contain loCorner?
          if (validBox.contains(loCorner))
            {
              found = true;
              const FArrayBox& thisLDF = levelLDF[dit];
              
              if (!m_useCubicInterp)
                {
                  Vector<Vector<Real> > linear(GLOBALDIM,Vector<Real>(2));
                  bool zeroWeight = false;
                  
                  for (int idir = 0; idir < GLOBALDIM; idir++)
                    {
                      Real f = fraction[idir];
                      
                      if (a_partialDerivative[idir] == 0)
                        {
                          linear[idir][0] = (-f + 1.0);
                          linear[idir][1] = ( f      );
                        }
                      else if (a_partialDerivative[idir] == 1)
                        {
                          linear[idir][0] = (-1.0);
                          linear[idir][1] = ( 1.0);
                        }
                      else
                        {
                          linear[idir][0] = (0.0);
                          linear[idir][1] = (0.0);
                          
                          zeroWeight = true;
                        }
                    }
                  
                  // If so, use multi-linear interpolation to get a return value
                  retval = 0.0;
                  
                  if (!zeroWeight)
                    {
                      // Iterate over all the corners
                      int comp = 0;
                      IVSIterator ivsit(m_cornersLinear);
                      for (ivsit.begin(); ivsit.ok(); ++ivsit)
                        {
                          IntVect curIV(loCorner);
                          IntVect incIV(ivsit());
                          Real curValue;
                          
                          // Index of the current data
                          curIV += incIV;
                          
                          // Check to see if current index in inside the data
                          if (dataBox.contains(curIV))
                            {
                              // If so, get the current data value
                              if (m_real_data != NULL)
                                {
                                  curValue = thisLDF(curIV,comp);
                                }
                              
                            }
                          else
                            {
                              // If not, use the no data value
                              curValue = m_noDataValue;
                            } 
                          
                          // Weight it appropriately based on the index fractions computed above
                          for (int idir = 0; idir < GLOBALDIM; idir++)
                            {
                              curValue *= linear[idir][incIV[idir]];
                            }
                          
                          // Add it into the total
                          retval += curValue;
                        }
                    }
                }
              else
                {
                  Vector<Vector<Real> > cubic(GLOBALDIM,Vector<Real>(4));
                  bool zeroWeight = false;
                  
                  for (int idir = 0; idir < GLOBALDIM; idir++)
                    {
                      Real f = fraction[idir];
                      
                      if (a_partialDerivative[idir] == 0)
                        {
                          cubic[idir][0] = (-0.5*f*f*f +     f*f - 0.5*f      );
                          cubic[idir][1] = ( 1.5*f*f*f - 2.5*f*f         + 1.0);
                          cubic[idir][2] = (-1.5*f*f*f + 2.0*f*f + 0.5*f      );
                          cubic[idir][3] = ( 0.5*f*f*f - 0.5*f*f              );
                        }
                      else if (a_partialDerivative[idir] == 1)
                        {
                          cubic[idir][0] = (-1.5*f*f + 2.0*f - 0.5);
                          cubic[idir][1] = ( 4.5*f*f - 5.0*f      );
                          cubic[idir][2] = (-4.5*f*f + 4.0*f + 0.5);
                          cubic[idir][3] = ( 1.5*f*f - 1.0*f      );
                        }
                      else if (a_partialDerivative[idir] == 2)
                        {
                          cubic[idir][0] = (-3.0*f + 2.0);
                          cubic[idir][1] = ( 9.0*f - 5.0);
                          cubic[idir][2] = (-9.0*f + 4.0);
                          cubic[idir][3] = ( 3.0*f - 1.0);
                        }
                      else if (a_partialDerivative[idir] == 3)
                        {
                          cubic[idir][0] = (-3.0);
                          cubic[idir][1] = ( 9.0);
                          cubic[idir][2] = (-9.0);
                          cubic[idir][3] = ( 3.0);
                        }
                      else
                        {
                          cubic[idir][0] = (0.0);
                          cubic[idir][1] = (0.0);
                          cubic[idir][2] = (0.0);
                          cubic[idir][3] = (0.0);
                          
                          zeroWeight = true;
                        }
                    }
                  
                  // If so, use multi-cubic interpolation to get a return value
                  retval = 0.0;
                  
                  if (!zeroWeight)
                    {
                      // Iterate over all the corners
                      int comp = 0;
                      IVSIterator ivsit(m_cornersCubic);
                      for (ivsit.begin(); ivsit.ok(); ++ivsit)
                        {
                          IntVect curIV(loCorner);
                          IntVect incIV(ivsit());
                          Real curValue;
                          
                          // Index of the current data
                          curIV += incIV;
                          
                          // Check to see if current index in insid the data
                          if (dataBox.contains(curIV))
                            {
                              // If so, get the current data value
                              curValue = thisLDF(curIV,comp);
                            }
                          else
                            {
                              // If not, use the no data value
                              curValue = m_noDataValue;
                            }
                          
                          // Weight it appropriate based on the index fractions computed above
                          for (int idir = 0; idir < GLOBALDIM; idir++)
                            {
                              curValue *= cubic[idir][incIV[idir]+1];
                            }
                          
                          // Add it into the total
                          retval += curValue;
                        }
                    }
                }
              
              if (a_partialDerivative.sum() == 0)
                {
                  // Adjust the level set from zero to m_value
                  retval -= m_value;
                }
              
              // Change the sign to change inside to outside
              if (!m_inside)
                {
                  retval = -retval;
                }
            } // end if this box contains the point
        } // end loop over grids on this level

      // if not found, drop down to the next level
      if (!found)
        {
          if (lev == 0)
            {
              // for now, Warn and set to default value
              MayDay::Warning("point not found in AMRDataFileIF::value()");
              retval = m_noDataValue;
            }
          else
            {
              lev -= 1;
              levelDomainBox.coarsen(m_refRatio[lev]);              
            }
        }
    }

  return retval;
}

BaseIF* AMRDataFileIF::newImplicitFunction() const
{
  AMRDataFileIF* dataFilePtr = new AMRDataFileIF(m_real_data,
                                                 m_char_data,
                                                 m_noDataValue,
                                                 m_baseDomain,
                                                 m_vectDx,
                                                 m_origin,
                                                 m_refRatio,
                                                 m_value,
                                                 m_inside,
                                                 m_useCubicInterp);

  return static_cast<BaseIF*>(dataFilePtr);
}


void AMRDataFileIF::ReadData(Real&             a_maxValue,
                             const char* const a_filename)
{
  // for now, assume origin is at zero
  m_origin = RealVect::Zero;

  // set up to read file
  string filename(a_filename);
  Vector<DisjointBoxLayout> vectGrids;
  Vector<LevelData<FArrayBox>* > vectData;
  Vector<string> vectNames;
  Real dx, dt, time;

  int status = ReadAMRHierarchyHDF5(filename,
                                    vectGrids,
                                    vectData,
                                    vectNames,
                                    m_baseDomain,
                                    dx,
                                    dt,
                                    time,
                                    m_refRatio,
                                    m_numLevels);

  CH_assert(status == 0);
  
  m_real_data = RefCountedPtr<Vector<LevelData<FArrayBox>* >>(new Vector<LevelData<FArrayBox>*>(m_numLevels, NULL));

  Box levelDomainBox;
  m_vectDx.resize(m_numLevels);
  
  for (int lev=0; lev<m_numLevels; lev++)
    {
      
      if (lev == 0)
        {
          m_vectDx[0] = dx*RealVect::Unit;
          levelDomainBox = m_baseDomain;
        }
      else
        {
          m_vectDx[lev] = m_vectDx[lev-1]/m_refRatio[lev-1];
          levelDomainBox = refine(levelDomainBox, m_refRatio[lev-1]);
        }
      
      // need ghost cells to handle coarse-fine BCs correctly
      IntVect ghostVect = IntVect::Unit;
      m_real_data->operator[](lev) = new LevelData<FArrayBox>(vectGrids[lev], vectNames.size(), ghostVect);
      // now copy
      vectData[lev]->copyTo(*(m_real_data->operator[](lev)));
      
      if (lev> 0)
        {
          // while we're here, might as well fill in ghost cells.
          //Use QuadCFInterp because it's more accurate
          QuadCFInterp interpolator(vectGrids[lev], &vectGrids[lev-1], m_vectDx[lev][0],
                                    m_refRatio[lev-1], vectNames.size(), levelDomainBox);

          interpolator.coarseFineInterp (*(m_real_data->operator[](lev)), *(m_real_data->operator[](lev-1)));
        }

      // basic extrapolation BCs for domain bc
      DataIterator dit = vectGrids[lev].dataIterator();
      
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box gridBox = vectGrids[lev][dit];
          FArrayBox& thisSDF = (*(m_real_data->operator[](lev)))[dit];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // low side
              Box edgeBoxLo = adjCellLo(gridBox, dir, 1);
              edgeBoxLo.grow(1);
              edgeBoxLo.grow(dir, -1);

              BoxIterator bitLo(edgeBoxLo);
              for (bitLo.begin(); bitLo.ok(); ++bitLo)
                {
                  IntVect iv = bitLo();
                  IntVect iv1 = iv + BASISV(dir);
                  IntVect iv2 = iv1 + BASISV(dir);
                  thisSDF(iv,0) = 2*thisSDF(iv1,0) - thisSDF(iv2,0);
                }

              // hi side
              Box edgeBoxHi = adjCellHi(gridBox, dir, 1);
              edgeBoxHi.grow(1);
              edgeBoxHi.grow(dir, -1);

              
              BoxIterator bitHi(edgeBoxHi);
              for (bitHi.begin(); bitHi.ok(); ++bitHi)
                {
                  IntVect iv = bitHi();
                  IntVect iv1 = iv - BASISV(dir);
                  IntVect iv2 = iv1 - BASISV(dir);
                  thisSDF(iv,0) = 2*thisSDF(iv1,0) - thisSDF(iv2,0);
                }
            } // end loop over directions
        } // end loop over grids
      
    } // end loop over levels

  // clean up
  for (int lev=0; lev<m_numLevels; lev++)
    {
      if (vectData[lev] != NULL)
        {
          delete (vectData[lev]);
          vectData[lev] = NULL;
        }
    }
  
}

void AMRDataFileIF::MakeCorners(void)
{
  // Make the end points of the linear interpolation box
  IntVect loLinear(0*IntVect::Unit);
  IntVect hiLinear(1*IntVect::Unit);

  // Make the unit cube box
  Box boxLinear(loLinear,hiLinear);

  // Save all its members (corners) as an IntVectSet
  m_cornersLinear.define(boxLinear);

  // Make the end points of the cubic interpolation box
  IntVect loCubic(-1*IntVect::Unit);
  IntVect hiCubic( 2*IntVect::Unit);

  // Make the unit cube box
  Box boxCubic(loCubic,hiCubic);

  // Save all its members (corners) as an IntVectSet
  m_cornersCubic.define(boxCubic);
}

#include "NamespaceFooter.H"
