!-------------------------------------------------------------------------
!> @brief     break a node and generate hierarchical nodes
!!
!> @param[in] Nod     - node number
!> @param[in] Kref    - refinement flag
!> @param[in] Iact    - T : generate active   sons
!!                      F : generate inactive sons
!> @date      Feb 2023
!-------------------------------------------------------------------------
!
subroutine nodbreak(Nod,Kref,Iact)
!
   use data_structure3D
!
   implicit none
!..Arguments
   integer, intent(in) :: Nod, Kref
   logical, intent(in) :: Iact
!..Local variables
   integer, dimension(27) :: ntype_sons, norder, nbcond, nsubd
   integer                :: nrsons, i, ison, icase
!
!-------------------------------------------------------------------------
!
#if HP3D_DEBUG
   integer :: iprint
   iprint = 0
   if (iprint.eq.1) then
      write(*,7001) Nod,Kref,Iact
 7001 format('nodbreak: Nod,Kref,Iact = ',i6,2x,i4,2x,l2)
   endif
#endif
!
!..record refinement kind
   NODES(Nod)%ref_kind=Kref
!
!..use Nod info to determine info about son nodes
   call set_break( NODES(Nod)%ntype,                       &
                   NODES(Nod)%ref_kind,                    &
                   NODES(Nod)%order,                       &
                   NODES(Nod)%bcond,                       &
                   NODES(Nod)%subd,                        &
                   nrsons, ntype_sons, norder, nbcond, nsubd )
!
!..generate the son nodes
   NODES(Nod)%nr_sons = nrsons
   icase = NODES(Nod)%case
!..Note: do not pass any member variable from NODES(Nod) into nodgen
!        if MAXNODS is increased, then NODES is reallocated
   do i=1,nrsons
      call nodgen( ntype_sons(i),                          &
                   icase,                                  &
                   nbcond(i),                              &
                   Nod,                                    &
                   norder(i),                              &
                   nsubd(i),                               &
                   Iact,                                   &
                   ison )
      if (i.eq.1) NODES(Nod)%first_son = ison
   enddo
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,7011) Nod
 7011 format('nodbreak: Nod ',i5,' HAS BEEN BROKEN')
      call pause
   endif
#endif
!
end subroutine nodbreak

