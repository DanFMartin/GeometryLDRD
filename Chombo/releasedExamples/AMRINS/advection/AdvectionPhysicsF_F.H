#ifndef _ADVECTIONPHYSICSF_F_H_
#define _ADVECTIONPHYSICSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRIEMANNF 
#define GUARDRIEMANNF 
// Prototype for Fortran procedure RIEMANNF ...
//
void FORTRAN_NAME( RIEMANNF ,riemannf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_FRA1(advVel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNF FORTRAN_NAME( inlineRIEMANNF, inlineRIEMANNF)
#define FORTNT_RIEMANNF FORTRAN_NAME( RIEMANNF, riemannf)

inline void FORTRAN_NAME(inlineRIEMANNF, inlineRIEMANNF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_FRA1(advVel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNF");
 FORTRAN_NAME( RIEMANNF ,riemannf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_FRA1(advVel)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNF 

#ifndef GUARDQUASILINEARUPDATE 
#define GUARDQUASILINEARUPDATE 
// Prototype for Fortran procedure QUASILINEARUPDATE ...
//
void FORTRAN_NAME( QUASILINEARUPDATE ,quasilinearupdate )(
      CHFp_FRA(AdWdx)
      ,CHFp_CONST_FRA(WHalf)
      ,CHFp_CONST_FRA1(cellVel)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_QUASILINEARUPDATE FORTRAN_NAME( inlineQUASILINEARUPDATE, inlineQUASILINEARUPDATE)
#define FORTNT_QUASILINEARUPDATE FORTRAN_NAME( QUASILINEARUPDATE, quasilinearupdate)

inline void FORTRAN_NAME(inlineQUASILINEARUPDATE, inlineQUASILINEARUPDATE)(
      CHFp_FRA(AdWdx)
      ,CHFp_CONST_FRA(WHalf)
      ,CHFp_CONST_FRA1(cellVel)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_QUASILINEARUPDATE");
 FORTRAN_NAME( QUASILINEARUPDATE ,quasilinearupdate )(
      CHFt_FRA(AdWdx)
      ,CHFt_CONST_FRA(WHalf)
      ,CHFt_CONST_FRA1(cellVel)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDQUASILINEARUPDATE 

}

#endif
