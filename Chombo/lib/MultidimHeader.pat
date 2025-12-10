#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM 1
#ifdef CH_USE1D 
#include "HEADERFILE"
#endif
#undef CH_SPACEDIM
#define CH_SPACEDIM 2
#ifdef CH_USE2D 
#include "HEADERFILE"
#endif
#undef CH_SPACEDIM
#define CH_SPACEDIM 3
#ifdef CH_USE3D 
#include "HEADERFILE"
#endif
#undef CH_SPACEDIM
#define CH_SPACEDIM 4
#ifdef CH_USE4D 
#include "HEADERFILE"
#endif
#undef CH_SPACEDIM
#define CH_SPACEDIM 5
#ifdef CH_USE5D 
#include "HEADERFILE"
#endif
#undef CH_SPACEDIM
#define CH_SPACEDIM 6
#ifdef CH_USE6D 
#include "HEADERFILE"
#endif
#undef CH_SPACEDIM
