#ifndef _ADVECTUTILF_F_H_
#define _ADVECTUTILF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDUDELS 
#define GUARDUDELS 
// Prototype for Fortran procedure UDELS ...
//
void FORTRAN_NAME( UDELS ,udels )(
      CHFp_FRA(uDotDelS)
      ,CHFp_FRA(sEdge)
      ,CHFp_FRA(vel)
      ,CHFp_BOX(gridBox)
      ,CHFp_REAL(dx)
      ,CHFp_INT(dir) );

#define FORT_UDELS FORTRAN_NAME( inlineUDELS, inlineUDELS)
#define FORTNT_UDELS FORTRAN_NAME( UDELS, udels)

inline void FORTRAN_NAME(inlineUDELS, inlineUDELS)(
      CHFp_FRA(uDotDelS)
      ,CHFp_FRA(sEdge)
      ,CHFp_FRA(vel)
      ,CHFp_BOX(gridBox)
      ,CHFp_REAL(dx)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_UDELS");
 FORTRAN_NAME( UDELS ,udels )(
      CHFt_FRA(uDotDelS)
      ,CHFt_FRA(sEdge)
      ,CHFt_FRA(vel)
      ,CHFt_BOX(gridBox)
      ,CHFt_REAL(dx)
      ,CHFt_INT(dir) );
}
#endif  // GUARDUDELS 

#ifndef GUARDADUALLSPEED 
#define GUARDADUALLSPEED 
// Prototype for Fortran procedure ADUALLSPEED ...
//
void FORTRAN_NAME( ADUALLSPEED ,aduallspeed )(
      CHFp_FRA(uDotDelS)
      ,CHFp_FRA(sEdge)
      ,CHFp_CONST_FRA(upEdge)
      ,CHFp_FRA(vel)
      ,CHFp_CONST_FRA(upCenter)
      ,CHFp_BOX(gridBox)
      ,CHFp_REAL(dx)
      ,CHFp_INT(dir) );

#define FORT_ADUALLSPEED FORTRAN_NAME( inlineADUALLSPEED, inlineADUALLSPEED)
#define FORTNT_ADUALLSPEED FORTRAN_NAME( ADUALLSPEED, aduallspeed)

inline void FORTRAN_NAME(inlineADUALLSPEED, inlineADUALLSPEED)(
      CHFp_FRA(uDotDelS)
      ,CHFp_FRA(sEdge)
      ,CHFp_CONST_FRA(upEdge)
      ,CHFp_FRA(vel)
      ,CHFp_CONST_FRA(upCenter)
      ,CHFp_BOX(gridBox)
      ,CHFp_REAL(dx)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_ADUALLSPEED");
 FORTRAN_NAME( ADUALLSPEED ,aduallspeed )(
      CHFt_FRA(uDotDelS)
      ,CHFt_FRA(sEdge)
      ,CHFt_CONST_FRA(upEdge)
      ,CHFt_FRA(vel)
      ,CHFt_CONST_FRA(upCenter)
      ,CHFt_BOX(gridBox)
      ,CHFt_REAL(dx)
      ,CHFt_INT(dir) );
}
#endif  // GUARDADUALLSPEED 

}

#endif
