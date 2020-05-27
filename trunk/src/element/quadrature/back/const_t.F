c----------------------------------------------------------------------
c
c   routine name       - const_T
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - routine initiates values of affine coordinates
c                        and weights for a Gaussian quadrature on
c                        a tetrahedron corresponding to Duffy's
c                        transformation
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine const_T
c
#include "syscom.blk"
#include "cint.blk"
c
c  ...force z to real type
      double precision z
c
      iprint=0
c
c  ...set up the number of integration points
c     defined in common as: NRGAUPO_TETS(7)
      do np=1,7
        NRGAUPO_TETS(np) = (np+1)*(np+1)*(np+2)
      enddo
c
c  ...set the integration rules for different orders p
c     defined in common as: NSELECT_TETS(7)
      do np=1,7
        NSELECT_TETS(np) = np
      enddo
c
c  ...set the quadrature data....
c
c  ...loop through quadrature rules
      do np=1,7
        if (iprint.eq.1) then
          write(*,7001) np
 7001     format('const_T: np = ',i2)
        endif
c
c  .....loop through integration points
        l=0
        ws = 0.d0
        do k=1,np+2
          z = XIGAUS1(k,np+2); wz = WAGAUS1(k,np+2)
          do j=1,np+1
            y = XIGAUS1(j,np+1); wy = WAGAUS1(j,np+1)
            do i=1,np+1
              x = XIGAUS1(i,np+1); wx = WAGAUS1(i,np+1)
              l=l+1
              XIGAUSS_TETS(1,l,np) = x*(1.d0-y)*(1.d0-z)
              XIGAUSS_TETS(2,l,np) = y*(1.d0-z)
              XIGAUSS_TETS(3,l,np) = z
              WAGAUSS_TETS(l,np) = wx*wy*wz*(1.d0-z)**2*(1.d0-y)
              ws = ws + WAGAUSS_TETS(l,np)
              if (iprint.eq.1) then
                write(*,7002) i,j,k,l,XIGAUSS_TETS(1:3,l,np)
 7002           format('const_T: i,j,k = ',3i3,' l = ',i4,
     .                 ' xi = ',3e12.5)
              endif
            enddo
          enddo
        enddo
        if (iprint.eq.1) then
          write(*,7003) ws
 7003     format('const_T: ws = ',e12.5)
          call pause
        endif
      enddo
c
c
ccc      return
c
c  ...verify the integration rule
c
c  ...loop through orders p
      do np=1,7
c
c  .....loop through monomials
        do i=0,2*np
        do j=0,2*np
        do k=0,2*np
          if (i+j+k.gt.2*np) cycle
c
c  .......evaluate the integral using 3D integration
          kint = NSELECT_TETS(np)
          nint= NRGAUPO_TETS(kint)
          val1 = 0.d0
          do l=1,nint
            x = XIGAUSS_TETS(1,l,kint)
            y = XIGAUSS_TETS(2,l,kint)
            z = XIGAUSS_TETS(3,l,kint)
            wa = WAGAUSS_TETS(l,kint)
            val1 = val1 + x**i*y**j*z**k*wa
          enddo
c
c  .......evaluate the integral using the 2D integration
          nord=(i+j+k+1)/2+1
          kint = NSELECT(nord)
          nint = NRGAUPO(kint)
          val2 = 0.d0
          do l=1,nint
            x = XIGAUSS(2,l,kint)
            y = XIGAUSS(3,l,kint)
            wa = WAGAUSS(l,kint)/2.d0
            val2 = val2 + x**i*y**j*(1.d0-x-y)**(k+1)/(k+1.d0)*wa
          enddo
          diff = abs(val1-val2)
          if (diff/abs(val1).gt.1.d-13) then
            write(*,7004) np
 7004       format('const_T: np = ',i2)
            write(*,7005) i,j,k
 7005       format('const_T: i,j,k = ',3i3)
            write(*,*) 'val1,val2 = ',val1,val2
            call pause
          endif
        enddo
        enddo
        enddo
      enddo
c
c
      end subroutine const_T
