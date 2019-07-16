!--------------------------------------------------------------------
!> Purpose : activate inactive nodes which are not constrained
!            but in the current mesh
!!
!> @date July 2019
!--------------------------------------------------------------------
!
subroutine refresh
!
   use data_structure3D
   use par_mesh  , only: DISTRIBUTED
   use MPI_param , only: RANK
!
   implicit none
!
   character(len=4) :: type
   integer, dimension(27) :: nodesl,norientl
   integer :: iprint,i,j,iel,nod,nfath,mdle,ibegin,iend, &
              nrsons,loc,subd
!
!--------------------------------------------------------------------
!
   iprint=0
!
   if (iprint.eq.1) then
      write(*,*) 'refresh: Begin'
   endif
!
!..reset visitation flags for all nodes
   call reset_visit
!  
!--------------------------------------------------------------------
! Step 1 : raise visitation flag for vertex, edge and face nodes of |
!          all active elements                                      |
!--------------------------------------------------------------------
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      call get_subd(mdle, subd)
      call elem_nodes(mdle, nodesl,norientl)
      type=NODES(mdle)%type
      ibegin=1
      iend  =nvert(type)+nedge(type)+nface(type)
      do i=ibegin,iend
         NODES(nodesl(i))%visit=1
!     ...if node is visited by an element within my subdomain,
!        add node to my subdomain (need its dofs). this flag will
!        indicate that dofs must be allocated in activation.
         if (DISTRIBUTED .and. (subd.eq.RANK)) then
            call set_subd(nodesl(i),subd)
         endif
      enddo
   enddo
!  
!--------------------------------------------------------------------
! Step 2: activate all inactive edge and face nodes whose father    |
!         node was not visited                                      |
!--------------------------------------------------------------------
!
!..loop over all nodes
   do nod=1,NRNODS
!
!  ...skip if a middle node
      select case(NODES(nod)%type)
         case('mdlb','mdln','mdlp','mdld') ; cycle
      end select
!
!  ...skip if the node has not been marked
      if (NODES(nod)%visit.eq.0) cycle
!
!  ...skip if active
      if (NODES(nod)%act.eq.1) cycle
!     
      if (iprint.eq.1) then
         write(*,7010) nod
 7010    format('refresh: INACTIVE MARKED NODE nod = ',i7)
      endif
!
      nfath=NODES(nod)%father
      if (nfath.le.0) then
         write(*,*) 'refresh: INCONSISTENCY: nod = ',nod
         stop
      endif
!
!  ...if father node has not been visited, activate the node
      if (NODES(nfath)%visit.eq.0) then
         call activate(nod)
         if (iprint.eq.1) then
            write(*,7020) nod
 7020       format('refresh: ACTIVATED nod = ',i6)
         endif
!
!     ...if this is the last son, deactivate the father
         nrsons = ubound(NODES(nfath)%sons,1)
         call locate(nod,NODES(nfath)%sons,nrsons, loc)
         if (loc.eq.nrsons) then
            call deactivate(nfath)
            if (iprint.eq.1) then
               write(*,7030) nfath
 7030          format('refresh: DEACTIVATED nfath, loc, nrsons = ',i6,2i3)
            endif
         endif
      endif
!
!..end of loop over nodes
   enddo
!
!..reset visitation flags
   call reset_visit
!
!
end subroutine refresh
