#include "typedefs.h"
subroutine cross_normal(Rn, Zcurl, Zcurj)
  implicit none
  real*8,dimension(3), intent(in)  :: Rn
  VTYPE, dimension(3), intent(in)  :: Zcurl
  VTYPE, dimension(3), intent(out) :: Zcurj
  !
  Zcurj(1) =   Rn(2)*Zcurl(3) - Rn(3)*Zcurl(2)
  Zcurj(2) =   Rn(3)*Zcurl(1) - Rn(1)*Zcurl(3)
  Zcurj(3) =   Rn(1)*Zcurl(2) - Rn(2)*Zcurl(1)
endsubroutine cross_normal
