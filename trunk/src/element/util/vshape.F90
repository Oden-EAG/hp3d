!----------------------------------------------------------------------
!> @brief computes vertex shape functions for a master triangle or quad
!!
!> @param[in]  Ftype - face type (TRIA,QUAD)
!> @param[in]  Xi    - master element coordinates
!> @param[out] Rlam  - vertex shape functions (i.e. baricentric coord.)
!> @param[out] Drlam - derivatives wrt master element coordinates
!!
!> @date Feb 2023
!----------------------------------------------------------------------
subroutine vshape2(Ftype,Xi, Rlam,Drlam)
!
   use parameters, only : MAXquadH !,MAXtriaH
   use node_types
!
   implicit none
!
   integer               , intent(in ) :: Ftype
   real(8),dimension(2  ), intent(in ) :: Xi
   real(8),dimension(  4), intent(out) :: Rlam
   real(8),dimension(2,4), intent(out) :: Drlam
!
   integer, dimension(5)          :: norder
   integer, dimension(4)          :: nedge_orient = (/0,0,0,0/)
   integer                        :: nrdof
   real(8), dimension(  MAXquadH) :: vshape
   real(8), dimension(2,MAXquadH) :: dvshape
!
   integer :: i,iprint
!----------------------------------------------------------------------
!
   iprint=0
!
!..compute vertex shape functions
   select case(Ftype)
      case(TRIA)
         norder=(/1,1,1,1,0/)
!         call shapeHt(Xi,norder(1:4),nedge_orient(1:3), nrdof, &
!                      vshape(1:MAXtriaH),dvshape(1:2,1:MAXtriaH))
      case(QUAD)
         norder=(/1,1,1,1,11/)
!         call shapeHq(Xi,norder,nedge_orient, nrdof,vshape,dvshape)
      case default
         write(*,7000) S_Type(Ftype)
 7000    format(' vshape2 : unknown figure! Ftype = ',a12)
   end select
   call shape2DH(Ftype,Xi,norder,nedge_orient, nrdof,vshape,dvshape)
!
!..store
   rlam(1:4)=vshape(1:4) ; drlam(1:2,1:4)=dvshape(1:2,1:4)
!
   if (iprint.eq.1) then
      write(*,7001) S_Type(Ftype)
 7001 format(' vshape2: Ftype = ',a4)
      do i=1,4
         write(*,7002) i,rlam(i),drlam(1:2,i)
 7002    format(' i = ',i1,'; rlam,drlam(1:2) = ',e12.5,4x,2(e12.5,2x))
      enddo
   endif
!
!
end subroutine vshape2
!
!----------------------------------------------------------------------
!> @brief computes vertex shape functions for a master block
!!
!> @param[in]  Etype - element type
!> @param[in]  Xi    - master element coordinates
!> @param[out] Rlam  - vertex shape functions (i.e. baricentric coord.)
!> @param[out] Drlam - derivatives wrt master element coordinates
!!
!> @date Feb 2023
!----------------------------------------------------------------------
subroutine vshape3(Etype,Xi, Rlam,Drlam)
!
   use parameters   , only : MAXbrickH
   use element_data , only : nedge
!
   implicit none
!
   integer                , intent(in ) :: Etype
   real(8), dimension(3  ), intent(in ) :: Xi
   real(8), dimension(  8), intent(out) :: Rlam
   real(8), dimension(3,8), intent(out) :: Drlam
!
   integer,dimension(19) :: norder
   integer,dimension(12) :: nedge_orient = (/0,0,0,0, 0,0,0,0, 0,0,0,0/)
   integer,dimension( 6) :: nface_orient = (/0,0,0,0, 0,0/)
   integer :: nrdof
   real(8) :: vshape(MAXbrickH)
   real(8) :: dvshape(3,MAXbrickH)
!
!----------------------------------------------------------------------
!
!..set up order of approximation
   call initiate_order(Etype, norder)
!
!..compute vertex shape functions
   call shape3DH(Etype,Xi,norder,nedge_orient,nface_orient, nrdof,vshape,dvshape)
!
!..store
   Rlam (    1:8) =  vshape(    1:8)
   Drlam(1:3,1:8) = dvshape(1:3,1:8)
!
end subroutine vshape3
