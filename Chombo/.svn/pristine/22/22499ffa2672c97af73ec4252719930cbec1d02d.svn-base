#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRDATAFILEIF_H_
#define _AMRDATAFILEIF_H_

#include <iostream>
#include <fstream>
using std::istream;
using std::ifstream;

#include "MayDay.H"
#include "RealVect.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "Vector.H"
#include "RefCountedPtr.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function reads data from a Chombo AMR HDF5 file and produces function values
    from the data read.  The data lies on an AMR hierarchy of rectilinear grids specified in the file.
    Function values within the data region are
    obtained using trilinear interpolation.  Function values outside the data
    region are the maximum value of the data read.

 */

class AMRDataFileIF: public BaseIF
{
public:
  ///

  ///
  /**
      Constructor specifying filename (a_filename), level set value (a_value), and whether
      inside the domain is where data is less than the level set value
      (a_inside).  Data is read from the file named, which must be a Chombo AMR HDF5 file.
   */
  AMRDataFileIF(const char* const           a_filename,
                const Real&                 a_value,
                const bool&                 a_inside,
                const bool&                 a_useCubicInterp = false);


  /// Copy constructor
  AMRDataFileIF(const AMRDataFileIF& a_inputIF);

  ///
  /**
      This is used by the factory (see below) to create a new object.  All
      objects created in this way share a refcounted pointer to their data.

      Constructor specifying a refcounted pointer to the data (a_real_data
      or a_char_data), the no data value (a_noDataValue), the spacing
      (a_spacing), the physical origin (a_origin), the number of grid points
      in each direction (a_num), a level set value (a_value), and whether
      inside the domain is where data is less than the level set value
      (a_inside).
   */
  AMRDataFileIF(const RefCountedPtr<Vector<LevelData<FArrayBox>* > > a_real_data,
                const RefCountedPtr<BaseFab<unsigned char> > a_char_data,
                const Real&                                  a_noDataValue,
                const Box&                                   a_baseDomain,
                const Vector<RealVect>&                      a_spacing,
                const RealVect&                              a_origin,
                const Vector<int>&                           a_refRatio,                                
                const Real&                                  a_value,
                const bool&                                  a_inside,
                const bool&                                  a_useCubicInterp = false);
  
  /// Destructor
  virtual ~AMRDataFileIF();

  ///
  /**
      Return the header information
   */
  virtual void GetHeader(Box&  a_baseDomain,
                         Vector<RealVect>& a_vectDx,
                         RealVect& a_origin) const;

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real& a_value,
                         bool& a_inside,
                         bool& a_useCubicInterp) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real& a_value,
                         const bool& a_inside,
                         const bool& a_useCubicInterp = false);

  ///
  /**
     value to use when we are outside
  */
  virtual void SetNoDataValue(const Real& a_value);

  ///
  /**
      Return the value of the function at a_point using trilinear interpolation
      of the data.  If a_point is outside the data then return the maximum data
      value.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  void GetFullHeader(Box&  a_baseDomain,
                     Vector<RealVect>& a_vectDx,
                     RealVect& a_origin)
  {
    a_baseDomain = m_baseDomain;
    a_vectDx = m_vectDx;
    a_origin = m_origin;
  }

  RefCountedPtr<FArrayBox> GetRealData(void)
  {
    return m_real_data;
  }

  RefCountedPtr<BaseFab<unsigned char> > GetCharData(void)
  {
    return m_char_data;
  }

protected:
  void ReadData(Real&             a_maxValue,
                const char* const a_filename);
  
  void MakeCorners(void);

  Box  m_baseDomain;       // number of grid points in each direction
  Vector<RealVect> m_vectDx;   // grid spacing (in physical coordinates)
  RealVect m_origin;    // grid origin (in physical coordinates)

  Real     m_value;     // level set value
  bool     m_inside;    // inside less than flag

  Vector<int> m_refRatio;
  int m_numLevels;
  
  // the data - copies all share this
  RefCountedPtr<Vector<LevelData<FArrayBox>* > >    m_real_data; // Real data stored as Real
  RefCountedPtr<BaseFab<unsigned char> > m_char_data; // Char data stored as unsigned char's

  Real m_noDataValue;   // no data value

  IntVectSet m_cornersLinear; // corners of a 1x1 box for linear interpolation
  IntVectSet m_cornersCubic;  // corners of a 4x4 box for cubic  interpolation

  bool m_useCubicInterp;

private:
  AMRDataFileIF()
  {
    MayDay::Abort("AMRDataFileIF uses strong construction");
  }

  void operator=(const AMRDataFileIF& a_inputIF)
  {
    MayDay::Abort("AMRDataFileIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
