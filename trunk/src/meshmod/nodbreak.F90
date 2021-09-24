!-------------------------------------------------------------------------
!> Purpose : break a node and generate hierarchical nodes
!!
!> @param[in] Nod     - node number
!> @param[in] Kref    - refinement flag
!> @param[in] Iact    - T : generate active   sons
!!                      F : generate inactive sons
!-------------------------------------------------------------------------
!
subroutine nodbreak(Nod,Kref,Iact)
!
   use data_structure3D
!
   implicit none
!..Arguments
   integer,               intent(in) :: Nod, Kref
   logical,               intent(in) :: Iact
!..Local variables
   character(len=4), dimension(27) :: type_sons
   integer,          dimension(27) :: norder, nbcond, nsubd
   integer                         :: nrsons, i, ison, icase
!
!-------------------------------------------------------------------------
!
#if DEBUG_MODE
   integer :: iprint
   iprint = 0
   if (iprint.eq.1) then
      write(*,7001) Nod,Kref,Iact
 7001 format('nodbreak: Nod,Kref,Iact = ',i6,2x,i4,2x,l2)
   endif
#endif
   select case(NODES(Nod)%type)
   case('medg','mdlt','mdlq','mdlb','mdlp','mdln','mdld')
   case default
     write(*,*) 'Nodbreak: Nod,NODES(Nod)%type = ',Nod,NODES(Nod)%type
     call result
   end select
!
!..record refinement kind
   NODES(Nod)%ref_kind=Kref
!
!..use Nod info to determine info about son nodes
   call set_break( NODES(Nod)%type,                        &
                   NODES(Nod)%ref_kind,                    &
                   NODES(Nod)%order,                       &
                   NODES(Nod)%bcond,                       &
                   NODES(Nod)%subd,                        &
                   nrsons, type_sons, norder, nbcond, nsubd )
!
!..generate the son nodes
   NODES(Nod)%nr_sons = nrsons
   icase = NODES(Nod)%case
!..Note: do not pass any member variable from NODES(Nod) into nodgen
!        if MAXNODS is increased, then NODES is reallocated
   do i=1,nrsons
      call nodgen( type_sons(i),                           &
                   icase,                                  &
                   nbcond(i),                              &
                   Nod,                                    &
                   norder(i),                              &
                   nsubd(i),                               &
                   Iact,                                   &
                   ison )
      if (i .eq. 1) NODES(Nod)%first_son = ison
   enddo
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7011) Nod
 7011 format('nodbreak: Nod ',i5,' HAS BEEN BROKEN')
      call pause
   endif
#endif
!
end subroutine nodbreak

