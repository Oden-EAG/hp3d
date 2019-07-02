!-------------------------------------------------------------------------
!> Purpose : element map of isoparametric element
!!
!! @param[in ] Mdle  - an (element) middle node number
!! @param[in ] Xi    - master element coordinates
!! @param[out] X     - physical coordinates
!! @param[in ] Dxdxi - derivative of physical coordinates
!!
!! @revision Nov 12
!-------------------------------------------------------------------------
!  Remark: the structure on an element routine is
!
!     determine integration points xiloc(:)
!     FOR i=1,#integration points
!       xi=xiloc(i)
!       compute geometry map at xi
!       ...
!     ENDLOOP
!
!  The first three instruction of this routine, marked with (*), do not 
!  depend upon Xi. Therefore, for a more efficient implementation of the 
!  element routine, they should be moved out of the loop. In other words,
!  rather then using this routine, it should be incorporated into the 
!  element routine, and a loop over integration points added where
!  indicated.
!-------------------------------------------------------------------------
!
subroutine geom_hp(Mdle,Xi, X,Dxdxi)
!
      use data_structure3D , only : NODES
      use parameters       , only : MAXbrickH
!
      implicit none
      integer,                 intent(in)  :: Mdle
      real*8,  dimension(3),   intent(in)  :: Xi
      real*8,  dimension(3),   intent(out) :: X
      real*8,  dimension(3,3), intent(out) :: Dxdxi
!
      integer,dimension(19) :: norder
      integer,dimension(12) :: nedge_orient
      integer,dimension( 6) :: nface_orient
      real*8,dimension(3,MAXbrickH) :: xnod
      real*8,dimension(  MAXbrickH) :: shapeH
      real*8,dimension(3,MAXbrickH) :: dshapeH
      integer :: k,i,nrdofH
!-------------------------------------------------------------------------
!
!  ...determine order of approximation, orientations, geometry dofs 
      call find_order( Mdle, norder)                                 !/
      call find_orient(Mdle, nedge_orient,nface_orient)              !/ (*)
      call nodcor(     Mdle, xnod)                                   !/
!
!////////  SET UP INTEGRATION POINTS  ///////////
!
!////////  LOOP OVER INTEGRATION POINTS  ////////
!
!  ...compute H1 shape functions
      call shape3H(NODES(Mdle)%Type,Xi,norder,nedge_orient,  &
                   nface_orient, nrdofH,shapeH,dshapeH)
!
!  ...geometry map 
      X(1:3)=0.d0 ; Dxdxi(1:3,1:3)=0.d0 
      do k=1,nrdofH
        X(1:3) = X(1:3) + xnod(1:3,k)*shapeH(k)
        do i=1,3
          Dxdxi(1:3,i) = Dxdxi(1:3,i) + xnod(1:3,k)*dshapeH(i,k)
        enddo
      enddo
! 
!////////  BODY OF ELEMENT ROUTINE //////////////
!
!////////  END OF LOOP  /////////////////////////
!      
!
endsubroutine geom_hp
