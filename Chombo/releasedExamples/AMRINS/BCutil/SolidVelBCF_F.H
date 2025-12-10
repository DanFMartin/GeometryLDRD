#ifndef _SOLIDVELBCF_F_H_
#define _SOLIDVELBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSOLIDVELBCF 
#define GUARDSOLIDVELBCF 
// Prototype for Fortran procedure SOLIDVELBCF ...
//
void FORTRAN_NAME( SOLIDVELBCF ,solidvelbcf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_REAL(bcVal)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_SOLIDVELBCF FORTRAN_NAME( inlineSOLIDVELBCF, inlineSOLIDVELBCF)
#define FORTNT_SOLIDVELBCF FORTRAN_NAME( SOLIDVELBCF, solidvelbcf)

inline void FORTRAN_NAME(inlineSOLIDVELBCF, inlineSOLIDVELBCF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_REAL(bcVal)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOLIDVELBCF");
 FORTRAN_NAME( SOLIDVELBCF ,solidvelbcf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_REAL(bcVal)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOLIDVELBCF 

}

#endif
