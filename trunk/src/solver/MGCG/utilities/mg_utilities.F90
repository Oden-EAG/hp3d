!  
!
#include "implicit_none.h"

   function l2norm(v)

   implicit none

   VTYPE,   intent(in) :: v(:)
   real*8 :: l2norm
#if C_MODE   
   real*8,  external   :: dznrm2
#else   
   real*8,  external   :: dnrm2
#endif   
   integer :: n
!
!------------------------------------------------------------
!
   n = size(v)

#if C_MODE
   l2norm = dznrm2(n,v,1)
#else
   l2norm = dnrm2(n,v,1)
#endif   

   end function l2norm
!
!------------------------------------------------------------
!------------------------------------------------------------
!
   function dotp(v,u)

   implicit none

   VTYPE, intent(in) :: v(:), u(:)
   VTYPE :: dotp
   integer :: n,m
#if C_MODE
   complex*16, external :: zdotc
#else
   real*8, external :: ddot
#endif      
!
!------------------------------------------------------------
!
   n = size(v); m = size(u)
   if (m .ne. n) then
      write(*,*) 'dotp: check sizes: size(v), size(u) = ', n, m
      stop 1
   endif

#if C_MODE
   dotp = zdotc(n,v,1,u,1)
#else
   dotp =  ddot(n,v,1,u,1)
#endif   

   end function dotp
