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
   use mpi_param , only: RANK
!
   implicit none
!
   character(4) :: type
   integer :: nodesl(27),norientl(27)
   integer :: mdlel(NRELES)
   integer :: iprint,i,j,iel,nod,nfath,mdle,ibegin,iend, &
              nrsons,loc,subd
!
!--------------------------------------------------------------------
!
   iprint=0
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,*) 'refresh: Begin'
   endif
#endif
!
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdlel(iel) = mdle
   enddo
!
!..reset visitation flags for all nodes
!$OMP PARALLEL
!$OMP DO
      do i=1,NRNODS
        NODES(i)%visit = 0
      enddo
!$OMP END DO
!  
!--------------------------------------------------------------------
! Step 1 : raise visitation flag for vertex, edge and face nodes of |
!          all active elements                                      |
!--------------------------------------------------------------------
!$OMP DO PRIVATE(mdle,subd,nodesl,norientl,type,ibegin,iend,i)
   do iel=1,NRELES
      mdle = mdlel(iel)
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
!$OMP END DO
!
!--------------------------------------------------------------------
! Step 2: activate all inactive edge and face nodes whose father    |
!         node was not visited                                      |
!--------------------------------------------------------------------
!
!..loop over all nodes
!$OMP DO SCHEDULE(DYNAMIC)       &
!$OMP PRIVATE(nfath,nrsons,loc)
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
#if DEBUG_MODE
!$OMP CRITICAL
      if (iprint.eq.1) then
         write(*,7010) nod
 7010    format('refresh: INACTIVE MARKED NODE nod = ',i7)
      endif
!$OMP END CRITICAL
#endif
!
      nfath=NODES(nod)%father
!
#if DEBUG_MODE
!$OMP CRITICAL
      if (nfath.le.0) then
         write(*,*) 'refresh: INCONSISTENCY: nod = ',nod
         stop
      endif
!$OMP END CRITICAL
#endif
!
!  ...if father node has not been visited, activate the node
      if (NODES(nfath)%visit.eq.0) then
         call activate(nod)
!
#if DEBUG_MODE
!$OMP CRITICAL
         if (iprint.eq.1) then
            write(*,7020) nod
 7020       format('refresh: ACTIVATED nod = ',i6)
         endif
!$OMP END CRITICAL
#endif
!
!     ...if this is the last son, deactivate the father
!         nrsons = ubound(NODES(nfath)%sons,1)
!         call locate(nod,NODES(nfath)%sons,nrsons, loc)
         nrsons = NODES(nfath)%nr_sons
         loc = nod - NODES(nfath)%first_son + 1
!         if (loc<0 .or. loc>nrsons) call pause
         if (loc.eq.nrsons) then
            call deactivate(nfath)
!
#if DEBUG_MODE
!$OMP CRITICAL
            if (iprint.eq.1) then
               write(*,7030) nfath
 7030          format('refresh: DEACTIVATED nfath, loc, nrsons = ',i6,2i3)
            endif
!$OMP END CRITICAL
#endif
!
         endif
      endif
!
!..end of loop over nodes
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
!
end subroutine refresh
