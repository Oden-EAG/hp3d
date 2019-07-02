!-------------------------------------------------------------------------
!> Purpose : break a node and generate hierarchical nodes
!!
!> @param[in] Nod     - node number
!> @param[in] Kref    - refinement flag
!> @param[in] Iact    - 1 : generate active   sons
!!                      0 : generate inactive sons
!> @param[in] Novert  - vertex nodes enclosing Nod (either edge, or face,
!!                                                  or element vertices) 
!> @param[in] Nr_vert - number of vertices
!-------------------------------------------------------------------------
!
subroutine nodbreak(Nod,Kref,Iact,Novert,Nr_vert)
!
  use data_structure3D
  implicit none
! Arguments
  integer,               intent(in) :: Nod, Kref, Iact, Nr_vert
  integer, dimension(8), intent(in) :: Novert
! Local variables
  real*8,           dimension(3)  :: x
  character(len=4), dimension(27) :: type
  integer,          dimension(27) :: norder, nfilter, nbcond 
  integer                         :: nrsons, iprint, i
!
!-------------------------------------------------------------------------
!
  iprint=0
!
! printing
  if (iprint.eq.1) then
     write(*,7001) Nod,Kref,Iact
7001 format('nodbreak: Nod,Kref,Iact = ',i6,2x,i4,2x,i2)
  endif
!
! record refinement kind
  NODES(Nod)%ref_kind=Kref 
!
! use the face or element centroid as an initial guess for the geoemtry dofs
  x(1:3)=0.d0
  if (Nr_vert.gt.0) then
     do i=1,Nr_vert
        x(1:3)=x(1:3)+NODES(Novert(i))%coord(1:3,1)
     enddo
     x(1:3)=x(1:3)/Nr_vert
  endif
!
! printing
  if (iprint.eq.1) then
     do i=1,Nr_vert
        write(*,7003) Novert(i), NODES(Novert(i))%coord(1:3,1)
     enddo
7003 format('nodbreak: nod, X = ',i6,3e12.5)
  endif
!
! use Nod info to determine info about son nodes
  nfilter(1:27) = 0
  call set_break( NODES(Nod)%type,                        &
                  NODES(Nod)%ref_kind,                    &
                  NODES(Nod)%ref_filter,                  &
                  NODES(Nod)%order,                       &
                  NODES(Nod)%bcond,                       &
                  nrsons, type, norder, nfilter, nbcond )
!
! generate the son nodes
  allocate(NODES(Nod)%sons(nrsons))
  do i=1,nrsons
     call nodgen( type(i),                                &
                  NODES(Nod)%case,                        &
                  nbcond(i),                              &
                  Nod,                                    &
                  norder(i),                              &
                  nfilter(i),                             &
                  Iact,                                   &
                  x,                                      &
                  NODES(Nod)%sons(i) )
  enddo
!
! printing
  if (iprint.eq.1) then
     write(*,7011) Nod 
7011 format('nodbreak: Nod ',i5,' HAS BEEN BROKEN')
     call pause
  endif
!
!
end subroutine nodbreak
