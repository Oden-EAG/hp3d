      real(8) function dreal_part(Z)
!
#include "syscom.blk"
!
#if C_MODE 
      dreal_part = real(Z)
#else
      dreal_part = Z
#endif
!
      end function dreal_part
!
!
!
      real(8) function dimag_part(Z)
!
#include "syscom.blk"
!
#if C_MODE
      dimag_part = aimag(Z)
#else
      dimag_part = 0.d0
#endif
!
end function dimag_part
