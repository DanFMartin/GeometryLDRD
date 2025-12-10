#ifndef _AVERAGEEDGEF_F_H_
#define _AVERAGEEDGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVERAGEEDGE 
#define GUARDAVERAGEEDGE 
// Prototype for Fortran procedure AVERAGEEDGE ...
//
void FORTRAN_NAME( AVERAGEEDGE ,averageedge )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_BOX(refBox) );

#define FORT_AVERAGEEDGE FORTRAN_NAME( inlineAVERAGEEDGE, inlineAVERAGEEDGE)
#define FORTNT_AVERAGEEDGE FORTRAN_NAME( AVERAGEEDGE, averageedge)

inline void FORTRAN_NAME(inlineAVERAGEEDGE, inlineAVERAGEEDGE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_BOX(refBox) )
{
 CH_TIMELEAF("FORT_AVERAGEEDGE");
 FORTRAN_NAME( AVERAGEEDGE ,averageedge )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(nRef)
      ,CHFt_BOX(refBox) );
}
#endif  // GUARDAVERAGEEDGE 

}

#endif
