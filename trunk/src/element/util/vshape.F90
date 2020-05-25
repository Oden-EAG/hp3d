!----------------------------------------------------------------------
!> Purpose : computes vertex shape functions for a master triangle or
!            quad      
!
!  @param[in ] Type  - shape      
!  @param[in ] Xi    - master element coordinates     
!  @param[out] Rlam  - vertex shape functions (i.e. baricentric coord.)
!  @param[out] Drlam - derivatives wrt master element coordinates
!----------------------------------------------------------------------
subroutine vshape2(Type,Xi, Rlam,Drlam)
!
   use parameters   , only : MAXquadH,MAXtriaH
!      
   implicit none
!
   character(len=4)      , intent(in ) :: Type
   real(8),dimension(2  ), intent(in ) :: Xi
   real(8),dimension(  4), intent(out) :: Rlam
   real(8),dimension(2,4), intent(out) :: Drlam
!
   integer,            dimension(5) :: norder
   integer, parameter, dimension(4) :: nedge_orient = (/0,0,0,0/)
   integer                          :: nrdof
   real(8), dimension(  MAXquadH)   :: vshape
   real(8), dimension(2,MAXquadH)   :: dvshape
!
   integer :: i,iprint
!----------------------------------------------------------------------
!      
   iprint=0
!      
!..compute vertex shape functions
   select case(Type)
      case('tria')
         norder=(/1,1,1,1,0/)
         call shapeHt(Xi,norder(1:4),nedge_orient(1:3), nrdof, &
                      vshape(1:MAXtriaH),dvshape(1:2,1:MAXtriaH))
      case('quad')
         norder=(/1,1,1,1,11/)
         call shapeHq(Xi,norder,nedge_orient, nrdof,vshape,dvshape)
      case default
         write(*,7000) Type
 7000    format(' vshape2 : unknown figure! Type = ',a12)
   end select
!      
!..store
   rlam(1:4)=vshape(1:4) ; drlam(1:2,1:4)=dvshape(1:2,1:4)
!
   if (iprint.eq.1) then
      write(*,7001) Type
 7001 format(' vshape2: Type = ',a4)
      do i=1,4
         write(*,7002) i,rlam(i),drlam(1:2,i)
 7002    format(' i = ',i1,'; rlam,drlam(1:2) = ',e12.5,4x,2(e12.5,2x))
      enddo
   endif
!
!
end subroutine vshape2
!
!
!----------------------------------------------------------------------
!> Purpose : computes vertex shape functions for a master block
!
!  @param[in ] Type  - block shape
!  @param[in ] Xi    - master element coordinates
!  @param[out] Rlam  - vertex shape functions (i.e. baricentric coord.)
!  @param[out] Drlam - derivatives wrt master element coordinates
!----------------------------------------------------------------------
subroutine vshape3(Type,Xi, Rlam,Drlam)
!
   use parameters   , only : MAXbrickH
   use element_data , only : nedge
!
   implicit none
!
   character(len=4)       , intent(in ) :: Type
   real(8), dimension(3  ), intent(in ) :: Xi
   real(8), dimension(  8), intent(out) :: Rlam
   real(8), dimension(3,8), intent(out) :: Drlam
!
   integer          ,dimension(19) :: norder
   integer,parameter,dimension(12) :: nedge_orient = (/0,0,0,0, 0,0,0,0, 0,0,0,0/)
   integer,parameter,dimension( 6) :: nface_orient = (/0,0,0,0, 0,0/)
   integer :: nrdof
   real(8) :: vshape(MAXbrickH)
   real(8) :: dvshape(3,MAXbrickH)
!
!----------------------------------------------------------------------
!
!..set up order of approximation
   call initiate_order(Type, norder)
!
!..compute vertex shape functions
   call shape3H(Type,Xi,norder,nedge_orient,nface_orient, nrdof,vshape,dvshape)
!
!..store
   rlam(1:8) = vshape(1:8); drlam(1:3,1:8) = dvshape(1:3,1:8)
!
!
end subroutine vshape3
!
