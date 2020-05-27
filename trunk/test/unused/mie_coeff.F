      program test_mie_coeff
      use parameters
#include "syscom.blk"
      complex*16 zref, z 
      complex*16 za(10), zb(10), zc(10), zd(10)
c
      write(*,*) 'test_mie:begin'
c
      n    = 5
      zref = cmplx(0.5, 0.0)
      z    = cmplx(0.5, 0.0)
      dmu  = 1.0
c
      call mie_coeff(n, zref, z, dmu, za,zb,zc,zd)
c
      do i=1,n
        write(*,7000)  i, za(i), zb(i), zc(i), zd(i)
 7000   format('i= ',i3, 4('(',e12.5,e12.5,')',2x))
      enddo
c
      write(*,*) 'test_mie:end'
c      
      end program test_mie_coeff
      
