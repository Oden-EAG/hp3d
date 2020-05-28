#include "typedefs.h"
subroutine compute_curl(Zdval, Zcurl)
  implicit none
  VTYPE, dimension(3,3), intent(in)  :: Zdval
  VTYPE, dimension(3),   intent(out) :: Zcurl
  !
  Zcurl(1) = Zdval(3,2) - Zdval(2,3)
  Zcurl(2) = Zdval(1,3) - Zdval(3,1)
  Zcurl(3) = Zdval(2,1) - Zdval(1,2)
endsubroutine compute_curl
