!----------------------------------------------------------------------
!
!   routine name       - curve_HermCur
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - Hermite parametrization curve
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
   subroutine curve_HermCur(No,Eta, X,Dxdeta)
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
      real(8) :: vshap(MAXP+1),dvshap(MAXP+1),ddvshap(MAXP+1)
!
      integer :: i,k,np,nord
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      if ((CURVES(No)%Type.ne.'HermCur').or.(NDIM.ne.3)) then
        write(*,7001)
 7001   format('curve_HermCur: WRONG INPUT')
        stop 1
      endif
!
!  ...get the curve endpoints
      do i=1,2
        np = CURVES(No)%EndPoNo(i)
        call pointr(np, gdof(1:3,i))
      enddo
!
!  ...get the tangent vectors and second derivatives
      gdof(1:3,3) = CURVES(No)%Rdata(1:3)
      gdof(1:3,4) = CURVES(No)%Rdata(4:6)
      gdof(1:3,5) = CURVES(No)%Rdata(7:9)
      gdof(1:3,6) = CURVES(No)%Rdata(10:12)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) No, Eta
 7002   format('curve_HermCur: No,Eta = ',i4,2x,f8.3)
        write(*,7003) (gdof(1:3,k),k=1,2)
        write(*,7003) (gdof(1:3,k),k=3,4)
        write(*,7003) (gdof(1:3,k),k=5,6)
 7003   format(6(3f8.3,2x))
      endif
#endif
!
!  ...determine C^2 polynomials
      nord=5
      call Gshape1(nord,Eta, vshap,dvshap,ddvshap)
!
      X(1:3) = 0.d0; Dxdeta(1:3) = 0.d0
      do k=1,6
        X(1:3) = X(1:3) + gdof(1:3,k)*vshap(k)
        Dxdeta(1:3) = Dxdeta(1:3) + gdof(1:3,k)*dvshap(k)
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7004)
 7004   format('curve_HermCur: X,Dxdeta = ')
        write(*,7003) X(1:3), Dxdeta(1:3)
        call pause
      endif
#endif
!
!
   end subroutine curve_HermCur
