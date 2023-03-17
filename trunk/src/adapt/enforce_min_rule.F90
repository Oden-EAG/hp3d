!----------------------------------------------------------------------
!
!   routine name       - enforce_min_rule
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine enforces the min rule for a FE mesh
!
!----------------------------------------------------------------------
!
   subroutine enforce_min_rule
!
   use data_structure3D
   use refinements
   use constrained_nodes
   use mpi_param,   only: ROOT, RANK
   use MPI,         only: MPI_COMM_WORLD, MPI_INTEGER, MPI_MAX, &
                          MPI_IN_PLACE
   use par_mesh,    only: DISTRIBUTED
   use bitvisit
   use par_ghost
!
   implicit none
!
!..work space for elem_nodes
   integer :: nodesl(27),norientl(27)
!
!..order for element nodes implied by the order of the middle node
   integer :: norder(19)
!
   integer :: ntype
!
   integer :: iel,mdle,nod,nrv,nre,nrf,nord
   integer :: i,j,k,nc,icase,nodp,nordh,nordv
   integer :: je,jf,ne1,ne2,ne3,ne4,is,nods
   integer :: nrsons,nordhs,nordvs
!
   integer, allocatable :: buffer(:)
!
   real(8) :: MPI_Wtime,start_time,end_time,t(2)
   integer :: ierr
!
#if DEBUG_MODE
   integer :: iprint
   iprint = 0
#endif
!
!----------------------------------------------------------------------
!
! TODO: this is only used here and in the MG solver at the minute;
!       may want to move into update_elem_order later though
!       (it is pretty inexpensive)
!..Get ghosted subdomain
   call get_ghost_subd
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..reset visitation flags
   call reset_visit
   call bitvisit_init(NRNODS)
!
!----------------------------------------------------------------------
!                 STEP 1: Minimum rule for faces
!----------------------------------------------------------------------
!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(iel,mdle,ntype,nrv,nre,nrf,nodesl,norientl,   &
!$OMP            norder,j,i,nod,nord,nc,icase,nodp)            &
!$OMP SCHEDULE(STATIC)
!..loop through elements in the current mesh
   do iel=1,NRELES_GHOST
      mdle = ELEM_GHOST(iel)
      ntype = NODES(Mdle)%ntype
      nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
      call get_connect_info(mdle, nodesl,norientl)
      call element_order(mdle,norientl, norder)
!
!  ...loop through edges and faces
      do j=1,nre+nrf
         i=nrv+j
         nod = nodesl(i); nord = norder(j)
!
         call save_min_order(mdle,nod,nord)
!
!     ...if a constrained node
         if (Is_inactive(nod)) then
!
!        ...identify the constraint case
            call decode2(NODES_CONSTR(i), nc,icase)
            select case(icase)
!
!        ...edge constrained by an edge
            case(11,12,37,38,47,48)
               nodp = NEDGC(nc)
               call save_min_order(mdle,nodp,nord)
!
!        ...horizontal edge constrained by a rectangular face
            case(26,28,33,36,63)
               nodp = NFACEC(nc)
               call save_min_order(mdle,nodp,nord*10+MAXP)
               do k=1,3,2
                  nodp = iabs(NFACE_CONS(k,nc))
                  call save_min_order(mdle,nodp,nord)
               enddo
!
!        ...vertical edge constrained by a rectangular face
            case(25,27,43,46,53)
               nodp = NFACEC(nc)
               call save_min_order(mdle,nodp,MAXP*10+nord)
               do k=2,4,2
                  nodp = iabs(NFACE_CONS(k,nc))
                  call save_min_order(mdle,nodp,nord)
               enddo
!
!        ...edge constrained by a triangular face
            case(75,76,77)
               nodp = NFACEC(nc)
               call save_min_order(mdle,nodp,nord)
               do k=1,3
                  nodp = iabs(NFACE_CONS(k,nc))
                  call save_min_order(mdle,nodp,nord)
               enddo
!
!        ...face node constrained by a face
            case(21,22,23,24,31,32,34,35,41,42,44,45,51,52,61,62,71,72,73,74)
               nodp = NFACEC(nc)
               call save_min_order(mdle,nodp,nord)
!
            case default
               !$OMP CRITICAL
               write(*,*) 'enforce_min_rule: icase = ',icase
               stop
               !$OMP END CRITICAL
            end select
         endif
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
!..Reduce to get nord of all nodes
!  (computing elem_nodes in first two stages is most expensive, this
!   way we only compute elem_nodes locally)
   if (DISTRIBUTED) then
      allocate(buffer(NRNODS)); buffer(:) = 0
!
!  ...loop through nodes touching my subdomain
      do nod=1,NRNODS
         if (visited(nod)) buffer(nod) = NODES(nod)%visit
      enddo
!
!  ...max used here to fill 0's; all procs should agree on shared nodes
      call MPI_Allreduce(MPI_IN_PLACE,buffer,NRNODS,MPI_INTEGER,MPI_MAX, MPI_COMM_WORLD,ierr)
!
      do nod=1,NRNODS
         if (visited(nod) .and. NODES(nod)%visit.ne.buffer(nod)) then
            write(*,*) 'enforce_min_rule: buffered nod order does not agree with original:', buffer(nod),' vs ', NODES(nod)%visit
            stop 2
         endif
         NODES(nod)%visit = buffer(nod)
      enddo
!
      deallocate(buffer)
   endif
!
   call bitvisit_finalize
!
!..loop through nodes
   do nod=1,NRNODS
!
!  ...active node
      if (Is_active(nod)) then
         nord = NODES(nod)%visit
         if (nord.eq.0) cycle
!
!     ...refined node, check order of its sons
         if (NODES(nod)%ref_kind.ne.0) then
            select case(NODES(nod)%ntype)
!
!        ...edge node
            case(MEDG)
               do is=1,2
                  nods = Son(nod,is)
                  nord = min(nord,NODES(nods)%visit)
               enddo
!
!        ...modify the node order
            call nodmod(nod, nord)
!
!        ...communicate the new order to the (inactive) sons
            do is=1,2
               nods = Son(nod,is)
               NODES(nods)%order = nord
            enddo
!
!        ...triangular face node
            case(MDLT)
               call nr_face_sons(MDLT,NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  nods = Son(nod,is)
                  select case(NODES(nods)%ntype)
                  case(MDLT)
                     nord = min(nord,NODES(nods)%visit)
                  case(MDLQ)
                     call decode(NODES(nods)%visit, nordh,nordv)
                     nord = min(nord,nordv)
                  end select
               enddo
!
!           ...modify the node order
               call nodmod(nod, nord)
!
!           ...communicate the new order to the (inactive) sons
               call nr_sons(MDLT,NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  nods = Son(nod,is)
                  select case(NODES(nods)%ntype)
                  case(MEDG); NODES(nods)%order = nord
                  case(MDLT); NODES(nods)%order = nord
                  case(MDLQ); NODES(nods)%order = nord*10+nord
                  end select
               enddo
!
!        ...rectangular face node
            case(MDLQ)
               call decode(nord, nordh,nordv)
               call nr_face_sons(MDLQ,NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  nods = Son(nod,is)
                  call decode(NODES(nods)%visit, nordhs,nordvs)
                  nordh = min(nordh,nordhs); nordv = min(nordv,nordvs)
               enddo
!
!           ...modify the node order
               nord = nordh*10+nordv
               call nodmod(nod, nord)
!
!           ...communicate the new order to the (inactive) mid-face sons
               do is=1,nrsons
                  !nods = NODES(nod)%sons(is)
                  nods = Son(nod,is)
                  NODES(nods)%order = nord
               enddo
!
!           ...communicate the new order to the (inactive) mid-edge sons
               select case(NODES(nod)%ref_kind)
               case(11)
                  do i=1,4
                     nods = Son(nod,nrsons+i)
                     select case(i)
                     case(1,3); NODES(nods)%order = nordv
                     case(2,4); NODES(nods)%order = nordh
                     end select
                  enddo
               case(10)
                  nods = Son(nod,nrsons+1)
                  NODES(nods)%order = nordv
               case(01)
                  nods = Son(nod,nrsons+1)
                  NODES(nods)%order = nordh
               end select
            end select
!
!     ...unrefined active node
         else
!
!        ...modify the node order
            nord = NODES(nod)%visit
            call nodmod(nod, nord)
!
!     ...if a refined node
         endif
!
!  ...if an active node
      endif
!
!..end of loop through nodes
   enddo
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (RANK .eq. ROOT) write(*,2020) end_time-start_time
2020 format(' enforce_min: ',f12.5,'  seconds')
!
!
contains
!
!
!> @date Mar 2023
!-----------------------------------------------------------------------
   subroutine save_min_order(Mdle,Nod,Nord)
!
      use data_structure3D
      use mpi_param
      implicit none
!
      integer, intent(in) :: Mdle,Nod,Nord
!
      integer :: nordh ,nordv
      integer :: nordh1,nordv1
      integer :: nordh2,nordv2
!
!  ...mark nodes connected to subdomain
      if (NODES(Mdle)%subd.eq.RANK) call visit(Nod)
!
!$OMP CRITICAL
      if (NODES(Nod)%visit.eq.0) then
         NODES(Nod)%visit = Nord
      else
         select case(NODES(Nod)%ntype)
         case(MEDG,MDLT)
            NODES(Nod)%visit = min(NODES(Nod)%visit,Nord)
         case(MDLQ)
            call decode(Nord, nordh1,nordv1)
            call decode(NODES(Nod)%visit, nordh2,nordv2)
            nordh = min(nordh1,nordh2)
            nordv = min(nordv1,nordv2)
            NODES(Nod)%visit = nordh*10+nordv
         end select
      endif
!$OMP END CRITICAL
!
   end subroutine save_min_order
!-----------------------------------------------------------------------
!
end subroutine enforce_min_rule
