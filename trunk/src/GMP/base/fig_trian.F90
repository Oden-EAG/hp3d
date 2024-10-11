!-----------------------------------------------------------------------
!> @brief Physical coordinates for a triangle parametrization, and
!!        their derivatives wrt to a GIVEN system of coordinates
!!
!> @param[in ] No      - triangle number
!> @param[in ] T       - coordinates of a point in the triangle
!!                       according to a GIVEN system
!> @param[in ] Norient - orientation of GLOBAL REFERENCE coordinates
!!                       wrt GIVEN coordinates
!> @param[out] X       - physical coordinates
!> @param[out] Dxdt    - derivatives of physical coordinates
!!
!> @date Nov 2012
!-----------------------------------------------------------------------
!
subroutine trian_local(No,T,Norient, X,Dxdt)
!
      use node_types, only: TRIA
      implicit none
      integer,               intent(in ) :: No,Norient
      real(8),dimension(2  ),intent(in ) :: T
      real(8),dimension(3  ),intent(out) :: X
      real(8),dimension(3,2),intent(out) :: Dxdt
!
      real(8),dimension(2  ) :: eta
      real(8),dimension(2,2) :: detadt
      real(8),dimension(3,2) :: dxdeta
      integer :: i
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!-----------------------------------------------------------------------
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7000) No,T,Norient
 7000   format(' trian_local: No,T,Norient = ',i6,2x,2(e12.5,2x,i1))
      endif
#endif
!
!  ...GIVEN -> GLOBAL REFERENCE
      call local2global(TRIA,T,Norient, eta,detadt)
!
!  ...GLOBAL REFERENCE -> PHYSICAL
      call trian(No,eta, X,dxdeta)
!
!  ...derivatives of PHYSICAL wrt GIVEN
      do i=1,2
        Dxdt(1:3,i) = dxdeta(1:3,1)*detadt(1,i) +  &
                      dxdeta(1:3,2)*detadt(2,i)
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        do i=1,3
          write(*,7001) i,X(i),Dxdt(i,1:2)
 7001     format(' i = ',i1,'; X(i),Dxdt(i,1:2) = ',3(e12.5,2x))
        enddo
      endif
#endif
!
end subroutine trian_local
!
!
!
!-----------------------------------------------------------------------
!> @brief Evaluates physical coordinates and its derivatives
!!        wrt to reference coordinates
!!
!> @param[in]  No     - a GMP triangle number
!> @param[in]  Eta    - reference coordinates of a point
!> @param[out] X      - physical coordinates of the point
!> @param[out] Dxdeta - derivatives of the physical coordinates
!!
!> @date Mar 2023
!-----------------------------------------------------------------------
!
subroutine trian(No,Eta, X,Dxdeta)
!
      use GMP , only : TRIANGLES , CURVES
!-----------------------------------------------------------------------
      implicit none
      integer,               intent(in)  :: No
      real(8),dimension(2),  intent(in)  :: Eta
      real(8),dimension(3),  intent(out) :: X
      real(8),dimension(3,2),intent(out) :: Dxdeta
!-----------------------------------------------------------------------
!  ...triangle vertices
      real(8),dimension(3,3) :: xv
!  ...edge curves
      integer,dimension(3) :: ncurv
      integer :: np, nc, i
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!-----------------------------------------------------------------------
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) No,Eta,TRIANGLES(No)%Type
 7001   format(' trian: No,Eta,Type = ',i8,2x,2f8.3,2x,a10)
      endif
#endif
!
      select case(TRIANGLES(No)%Type)
!
!  ...linear or plane triangle
      case('PlaneTri')
!  .....get the edges and then the vertices
        do i=1,3
          ncurv(i)=TRIANGLES(No)%EdgeNo(i)
          np=TRIANGLES(No)%VertNo(i) ; call pointr(np, xv(1:3,i))
        enddo
!  .....check compatibility
        do i=1,3
          nc=abs(ncurv(i))
          if (CURVES(nc)%Type.ne.'Seglin') then
            write(*,7002) No,i,nc,CURVES(nc)%Type
 7002       format(' trian: INCOMPATIBLE DEFINITIONS: No,i,ic=',  &
                    i4,2x,i2,i4,' type = ',a10)
            stop
          endif
        enddo
!  .....accumulate
        do i=1,3
          Dxdeta(i,1) = xv(i,2) - xv(i,1)
          Dxdeta(i,2) = xv(i,3) - xv(i,1)
          X(i) = xv(i,1) + Eta(1)*Dxdeta(i,1) + Eta(2)*Dxdeta(i,2)
        enddo
!
!  ...transfinite interpolation triangle................................
      case('TransTri') ; call trian_TransTri(No,Eta, X,Dxdeta)
!
!  ...parametric transfinite interpolation triangle.....................
      case('PTITri')   ; call trian_PTITri(  No,Eta, X,Dxdeta)
!
!  ...implicit triangle (CURRENTLY NOT WORKING!)........................
      !case('ImpliTri') ; call trian_ImpliTri(No,Eta, X,Dxdeta)
!
!  ...spherical octant (LEGACY).........................................
      !case('SpherTri') ; call trian_SpherTri(No,Eta, X,Dxdeta)
!
!  ...quarter of a circle (LEGACY)......................................
      !case('QtCirTri') ; call trian_QtCirTri(No,Eta, X,Dxdeta)
!
!  ...a part of spherical octant (LEGACY)...............................
      !case('PaSphTri') ; call trian_PaSphTri(No,Eta, X,Dxdeta)
!
!  ...G1 reconstructed surface triangle (LEGACY)........................
      !case('G1RecTri') ; call trian_G1RecTri(No,Eta, X,Dxdeta)
!
!  ...cylindrical coordinates triangle
      case('CylTri')   ; call trian_CylTri(No,Eta, X,Dxdeta)
!
      case default
      write(*,7003) TRIANGLES(No)%Type
 7003 format(' trian: unknown Type = ',a10)
      stop
!
      end select
!
!
end subroutine trian
