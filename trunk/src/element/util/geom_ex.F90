subroutine geom_ex(Mdle,Xi, X,Dxdxi)
!
      integer,                  intent(in)  :: Mdle
      real(8),  dimension(3),   intent(in)  :: Xi
      real(8),  dimension(3),   intent(out) :: X
      real(8),  dimension(3,3), intent(out) :: Dxdxi
!
!  ...redirect to old routine to ensure backward compatibility
      call exact_geom(Mdle,Xi, X,Dxdxi)
!
end subroutine geom_ex
!
!
!------------------------------------------------------------------------
!> @brief routine computes physical coordinates for a point in the
!!           master element using the exact geometry maps
!!
!> @param[in]  Mdle  - middle node
!> @param[in]  Xi    - a point within the corresponding master element
!> @param[out] X     - coordinates of the corresponding physical point
!> @param[out] Dxdxi - derivatives wrt master element coordinates
!!
!> @date Feb 2023
!------------------------------------------------------------------------
subroutine exact_geom(Mdle,Xi, X,Dxdxi)
!
      use data_structure3D
!
      implicit none
      integer,                 intent(in)  :: Mdle
      real(8), dimension(3),   intent(in)  :: Xi
      real(8), dimension(3),   intent(out) :: X
      real(8), dimension(3,3), intent(out) :: Dxdxi
!
      real(8), dimension(  8) ::  shapH
      real(8), dimension(3,8) :: dshapH

      real(8) :: eta(3),dxdeta(3,3),detadxi(3,3),etav(3,8)
      integer :: iflag, i, j, k, no
      integer :: ntype
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!------------------------------------------------------------------------
!
      ntype=NODES(Mdle)%ntype
!
!  ...element vertices in the reference space
      call refel(Mdle, iflag,no,etav)
!
!  ...vertex shape functions
      call vshape3(ntype,Xi, shapH,dshapH)
!
!  ...determine refinement map : Eta = Eta(Xi)
      eta(1:3)=0.d0 ; detadxi(1:3,1:3)=0.d0
      do k=1,nvert(ntype)
        eta(1:3) = eta(1:3) + etav(1:3,k)*shapH(k)
        do i=1,3
          detadxi(1:3,i) = detadxi(1:3,i) + etav(1:3,k)*dshapH(i,k)
        enddo
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Mdle,iflag,no
7001    format('exact_geom: Mdle = ',i5,' iflag,no  = ',2i4,' etav = ')
        do i=1,8
          write(*,7002) i,etav(1:3,i)
7002      format(i3,2x,3f8.3)
        enddo
        write(*,7003) eta(1:3)
7003    format('            eta = ',3e12.5)
        call pause
      endif
#endif
!
!  ...compose refinement map with GMP map
      select case(iflag)
        case(5) ; call prism(no,eta, x,dxdeta)
        case(6) ; call  hexa(no,eta, x,dxdeta)
        case(7) ; call tetra(no,eta, x,dxdeta)
        case(8) ; call pyram(no,eta, x,dxdeta)
        case default
          write(*,*) 'exact_geom: Mdle, type, iflag = ',Mdle,S_Type(ntype),iflag
          call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
      endselect
!
!  ...adjust derivatives
      Dxdxi(1:3,1:3)=0.d0
      do i=1,3
        do j=1,3
          Dxdxi(1:3,i) = Dxdxi(1:3,i) + dxdeta(1:3,j)*detadxi(j,i)
        enddo
      enddo
!
end subroutine exact_geom
