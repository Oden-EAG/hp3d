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
!
   implicit none
!..Arguments
   integer,               intent(in) :: Nod, Kref, Iact, Nr_vert
   integer, dimension(8), intent(in) :: Novert
!..Local variables
   character(len=4), dimension(27) :: type_sons
   integer,          dimension(27) :: norder, nfilter, nbcond, nsubd
   integer                         :: nrsons, iprint, i, ison
!
!-------------------------------------------------------------------------
!
   iprint=0
!
!..printing
   if (iprint.eq.1) then
      write(*,7001) Nod,Kref,Iact
 7001 format('nodbreak: Nod,Kref,Iact = ',i6,2x,i4,2x,i2)
   endif
!
!..record refinement kind
   NODES(Nod)%ref_kind=Kref
!
!..use Nod info to determine info about son nodes
   nfilter(1:27) = 0
   call set_break( NODES(Nod)%type,                        &
                   NODES(Nod)%ref_kind,                    &
                   NODES(Nod)%ref_filter,                  &
                   NODES(Nod)%order,                       &
                   NODES(Nod)%bcond,                       &
                   NODES(Nod)%subd,                        &
                   nrsons, type_sons, norder, nfilter, nbcond, nsubd )
!
!..generate the son nodes
!   allocate(NODES(Nod)%sons(nrsons))
   NODES(Nod)%nr_sons = nrsons
   do i=1,nrsons
      call nodgen( type_sons(i),                           &
                   NODES(Nod)%case,                        &
                   nbcond(i),                              &
                   Nod,                                    &
                   norder(i),                              &
                   nfilter(i),                             &
                   nsubd(i),                               &
                   Iact,                                   &
                   ison )
!                   NODES(Nod)%sons(i) )
      if (i .eq. 1) NODES(Nod)%first_son = ison
   enddo
!
!..printing
   if (iprint.eq.1) then
      write(*,7011) Nod
 7011 format('nodbreak: Nod ',i5,' HAS BEEN BROKEN')
      call pause
   endif
!
!
end subroutine nodbreak
