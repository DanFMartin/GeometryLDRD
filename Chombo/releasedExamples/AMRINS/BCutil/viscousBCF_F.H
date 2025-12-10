#ifndef _VISCOUSBCF_F_H_
#define _VISCOUSBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDVISCOUSBC 
#define GUARDVISCOUSBC 
// Prototype for Fortran procedure VISCOUSBC ...
//
void FORTRAN_NAME( VISCOUSBC ,viscousbc )(
      CHFp_FRA(viscous)
      ,CHFp_BOX(gridBox)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(hiBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(nComp) );

#define FORT_VISCOUSBC FORTRAN_NAME( inlineVISCOUSBC, inlineVISCOUSBC)
#define FORTNT_VISCOUSBC FORTRAN_NAME( VISCOUSBC, viscousbc)

inline void FORTRAN_NAME(inlineVISCOUSBC, inlineVISCOUSBC)(
      CHFp_FRA(viscous)
      ,CHFp_BOX(gridBox)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(hiBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(nComp) )
{
 CH_TIMELEAF("FORT_VISCOUSBC");
 FORTRAN_NAME( VISCOUSBC ,viscousbc )(
      CHFt_FRA(viscous)
      ,CHFt_BOX(gridBox)
      ,CHFt_BOX(loBox)
      ,CHFt_BOX(hiBox)
      ,CHFt_INT(dir)
      ,CHFt_INT(nComp) );
}
#endif  // GUARDVISCOUSBC 

}

#endif
