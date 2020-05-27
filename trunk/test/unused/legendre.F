      program test_legendre

#include "cint.blk"
      
      real*8 :: xi, w, err_val, product
      double precision, dimension(0:20) :: poly, dpoly
      double precision, dimension(1:20) :: val
      integer :: i, j, k, nint, np

      poly =0.d0; dpoly= 0.d0;
      call legendre_I(  0.3d0, 6, poly, dpoly )
      write(*,*) ' 0.3 poly = ', poly(1:8)
      write(*,*) ' 0.3dpoly = ', dpoly(1:8)      

      stop 1
c     Test up to 6
      np = 6
      
c     combination of legendre polynomial
      do i=1, np
        do j=1, np

c     ....number of integration points
          nint = i+j + 1 
c     
          product = 0.0
c     ....loop through integration points
          do k=1, nint
c     ......rescaling
            xi = 2.0*XIGAUS1(k, nint)-1.0
            w  = 2.0*WAGAUS1(k, nint)
            
            call legendre( xi, 6, poly, dpoly )
            product = product + w*poly(i)*poly(j)
          enddo
c     ....orthogonality condition
          if (i.eq.j) then
            err_val = abs(product -2.0/(2.0*i+1.0))
          else
            err_val = abs(product)
          endif
c          
          write(*,*) 'i=',i,'j=',j,'err=',err_val
        enddo
      enddo
      
      end program test_legendre
      
