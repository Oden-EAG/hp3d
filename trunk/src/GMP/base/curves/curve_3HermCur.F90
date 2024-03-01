!----------------------------------------------------------------------
!
!   routine name       - curve_3HermCur
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - cubic Hermite curve
!
!   arguments :
!     in:
!               No     - the curve number
!               Eta    - reference coordinate  (between 0 and 1)
!     out:
!               X      - physical coordinates of the point
!               Dxdeta - derivatives of the physical coordinates wrt
!                        reference coordinate
!
!----------------------------------------------------------------------
!
   subroutine curve_3HermCur(No,Eta, X,Dxdeta)
!
      use GMP
      use parameters , only : MAXP
!
      implicit none
!
      integer :: No
      real(8) :: Eta,X(3),Dxdeta(3)
!
!  ...curve dof
      real(8) :: gdof(1:3,6)
!
!  ...fifth order C^2-conforming shape functions
      real(8) :: vshap(MAXP+1),dvshap(MAXP+1)
!
      integer :: i,k,np
!
      integer :: iprint
      iprint=0
!
      if ((CURVES(No)%Type.ne.'3HermCur').or.(NDIM.ne.3)) then
        write(*,7001)
 7001   format('curve_3HermCur: WRONG INPUT')
        stop
      endif
!
!  ...get curve endpoints
      do i=1,2
        np = CURVES(No)%EndPoNo(i)
        call pointr(np, gdof(1:3,i))
      enddo
!
!  ...get velocities at endpoints
      gdof(1:3,3) = CURVES(No)%Rdata(1:3)
      gdof(1:3,4) = CURVES(No)%Rdata(4:6)
!!!      gdof(1:3,5) = CURVES(No)%Rdata(7:9)
!!!      gdof(1:3,6) = CURVES(No)%Rdata(10:12)
      if (iprint.eq.1) then
        write(*,7002) No, Eta
 7002   format('curve_HermCur: No,Eta = ',i4,2x,f8.3)
        write(*,7003) (gdof(1:3,k),k=1,2)
        write(*,7003) (gdof(1:3,k),k=3,4)
!!!        write(*,7003) (gdof(1:3,k),k=5,6)
 7003   format(4(3f8.3,2x))
      endif
!
!  ...determine polynomials
      call cubic_Hermite(Eta, vshap,dvshap)
!
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
!
!
   end subroutine curve_3HermCur
!
!
!----------------------------------------------------------------------
!
   subroutine cubic_Hermite(Eta, Vshap,Dvshap)
!
!----------------------------------------------------------------------
!
      implicit none
!
      real(8) :: Eta
      real(8) :: Vshap(4),Dvshap(4)
!
      Vshap(1) = 2.d0*Eta**3 - 3.d0*Eta**2 + 1.d0
      Vshap(2) = -2.d0*Eta**3 + 3.d0*Eta**2
      Vshap(3) = Eta**3 - 2.d0*Eta**2 + Eta
      Vshap(4) = Eta**3 - Eta**2
!
      Dvshap(1) = 6.d0*Eta**2 - 6.d0*Eta
      Dvshap(2) = -6.d0*Eta**2 + 6.d0*Eta
      Dvshap(3) = 3.d0*Eta**2 - 4.d0*Eta**1 + 1.d0
      Dvshap(4) = 3.d0*Eta**2 - 2.d0*Eta
!
!
   end subroutine cubic_Hermite
