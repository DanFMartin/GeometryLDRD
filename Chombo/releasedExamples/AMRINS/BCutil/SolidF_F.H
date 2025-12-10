#ifndef _SOLIDF_F_H_
#define _SOLIDF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSOLIDBCF 
#define GUARDSOLIDBCF 
// Prototype for Fortran procedure SOLIDBCF ...
//
void FORTRAN_NAME( SOLIDBCF ,solidbcf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_REAL(bcVal)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_SOLIDBCF FORTRAN_NAME( inlineSOLIDBCF, inlineSOLIDBCF)
#define FORTNT_SOLIDBCF FORTRAN_NAME( SOLIDBCF, solidbcf)

inline void FORTRAN_NAME(inlineSOLIDBCF, inlineSOLIDBCF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_REAL(bcVal)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOLIDBCF");
 FORTRAN_NAME( SOLIDBCF ,solidbcf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_REAL(bcVal)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOLIDBCF 

#ifndef GUARDSOLIDEXTRAPBCF 
#define GUARDSOLIDEXTRAPBCF 
// Prototype for Fortran procedure SOLIDEXTRAPBCF ...
//
void FORTRAN_NAME( SOLIDEXTRAPBCF ,solidextrapbcf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_REAL(bcVal)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_SOLIDEXTRAPBCF FORTRAN_NAME( inlineSOLIDEXTRAPBCF, inlineSOLIDEXTRAPBCF)
#define FORTNT_SOLIDEXTRAPBCF FORTRAN_NAME( SOLIDEXTRAPBCF, solidextrapbcf)

inline void FORTRAN_NAME(inlineSOLIDEXTRAPBCF, inlineSOLIDEXTRAPBCF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_REAL(bcVal)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOLIDEXTRAPBCF");
 FORTRAN_NAME( SOLIDEXTRAPBCF ,solidextrapbcf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_REAL(bcVal)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOLIDEXTRAPBCF 

#ifndef GUARDSLOPEBCSF 
#define GUARDSLOPEBCSF 
// Prototype for Fortran procedure SLOPEBCSF ...
//
void FORTRAN_NAME( SLOPEBCSF ,slopebcsf )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(loVal)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_CONST_REAL(hiVal)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi) );

#define FORT_SLOPEBCSF FORTRAN_NAME( inlineSLOPEBCSF, inlineSLOPEBCSF)
#define FORTNT_SLOPEBCSF FORTRAN_NAME( SLOPEBCSF, slopebcsf)

inline void FORTRAN_NAME(inlineSLOPEBCSF, inlineSLOPEBCSF)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(loVal)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_CONST_REAL(hiVal)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi) )
{
 CH_TIMELEAF("FORT_SLOPEBCSF");
 FORTRAN_NAME( SLOPEBCSF ,slopebcsf )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(loVal)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_CONST_REAL(hiVal)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi) );
}
#endif  // GUARDSLOPEBCSF 

}

#endif
