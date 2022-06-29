c----------------------------------------------------------------------
c
c   routine name       - curve_3HermCur
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - cubic Hermite curve
c
c   arguments :
c     in:
c               No     - the curve number
c               Eta    - reference coordinate  (between 0 and 1)
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        reference coordinate
c
c----------------------------------------------------------------------
c
      subroutine curve_3HermCur(No,Eta, X,Dxdeta)
c
      use GMP
      use parameters , only : MAXP
#include "syscom.blk"
c
      dimension X(3),Dxdeta(3)
c
c  ...curve dof
      dimension gdof(1:3,6)
c
c  ...fifth order C^2-conforming shape functions
      dimension vshap(MAXP+1),dvshap(MAXP+1)
c
      iprint=0
c
      if ((CURVES(No)%Type.ne.'3HermCur').or.(NDIM.ne.3)) then
        write(*,7001)
 7001   format('curve_3HermCur: WRONG INPUT')
        stop
      endif
c
c  ...get curve endpoints
      do i=1,2
        np = CURVES(No)%EndPoNo(i)
        call pointr(np, gdof(1:3,i))
      enddo
c
c  ...get velocities at endpoints
      gdof(1:3,3) = CURVES(No)%Rdata(1:3)
      gdof(1:3,4) = CURVES(No)%Rdata(4:6)
ccc      gdof(1:3,5) = CURVES(No)%Rdata(7:9)
ccc      gdof(1:3,6) = CURVES(No)%Rdata(10:12)
      if (iprint.eq.1) then
        write(*,7002) No, Eta
 7002   format('curve_HermCur: No,Eta = ',i4,2x,f8.3)
        write(*,7003) (gdof(1:3,k),k=1,2)
        write(*,7003) (gdof(1:3,k),k=3,4)
ccc        write(*,7003) (gdof(1:3,k),k=5,6)
 7003   format(4(3f8.3,2x))
      endif
c
c  ...determine polynomials
      call cubic_Hermite(Eta, vshap,dvshap)
c
      X(1:3) = 0.d0; Dxdeta(1:3) = 0.d0
      do k=1,4
        X(1:3) = X(1:3) + gdof(1:3,k)*vshap(k)
        Dxdeta(1:3) = Dxdeta(1:3) + gdof(1:3,k)*dvshap(k)
      enddo
      if (iprint.eq.1) then
        write(*,7004)
 7004   format('curve_3HermCur: X,Dxdeta = ')
        write(*,7003) X(1:3), Dxdeta(1:3)
        call pause
      endif
c
c
      end
c
c
c----------------------------------------------------------------------
c
      subroutine cubic_Hermite(Eta, Vshap,Dvshap)
c
c----------------------------------------------------------------------
c
#include "syscom.blk"
c
      dimension Vshap(4),Dvshap(4)
c
      Vshap(1) = 2.d0*Eta**3 - 3.d0*Eta**2 + 1.d0
      Vshap(2) = -2.d0*Eta**3 + 3.d0*Eta**2
      Vshap(3) = Eta**3 - 2.d0*Eta**2 + Eta
      Vshap(4) = Eta**3 - Eta**2
c
      Dvshap(1) = 6.d0*Eta**2 - 6.d0*Eta
      Dvshap(2) = -6.d0*Eta**2 + 6.d0*Eta
      Dvshap(3) = 3.d0*Eta**2 - 4.d0*Eta**1 + 1.d0
      Dvshap(4) = 3.d0*Eta**2 - 2.d0*Eta
c
c
      end
