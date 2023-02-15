!----------------------------------------------------------------------
!
!   routine name       - enforce_min_rule
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 23
!
!   purpose            - routine enforces the min rule for a FE mesh
!
!
!----------------------------------------------------------------------
!
   subroutine enforce_min_rule
!
   use data_structure3D
   use refinements
   use constrained_nodes
   use mpi_param  , only: ROOT,RANK
   use MPI        , only: MPI_COMM_WORLD,MPI_COMM_WORLD
!
   implicit none
!
!..work space for elem_nodes
   integer :: nodesl(27),norientl(27)
!
!..order for element nodes implied by the order of the middle node
   integer :: norder(19)
!
   character(len=4) :: etype
!
   integer :: iel,mdle,nod,nrv,nre,nrf,nord
   integer :: i,j,k,nc,icase,nodp,nordh,nordv
   integer :: je,jf,ne1,ne2,ne3,ne4,is,nods
   integer :: nrsons,nordhs,nordvs
!
   real(8) :: MPI_Wtime,start_time,end_time
   integer :: ierr
!
#if DEBUG_MODE
   integer :: iprint
   iprint = 0
#endif
!
!----------------------------------------------------------------------
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..reset visitation flags
   call reset_visit
!
!----------------------------------------------------------------------
!  STEP 1:  Save order implied by middle nodes to regular nodes
!           and parents of constrained nodes
!----------------------------------------------------------------------
!
!..loop through elements in the current mesh
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
      call get_connect_info(mdle, nodesl,norientl)
      call element_order(mdle,norientl, norder)
!
!  ...loop through edges and faces
      do j=1,nre+nrf
         i=nrv+j
         nod = nodesl(i); nord = norder(j)
!
         call save_min_order(nod,nord)
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
               call save_min_order(nodp,nord)
!
!        ...horizontal edge constrained by a rectangular face
            case(26,28,33,36,63)
               nodp = NFACEC(nc)
               call save_min_order(nodp,nord*10+MAXP)
               do k=1,3,2
                  nodp = iabs(NFACE_CONS(k,nc))
                  call save_min_order(nodp,nord)
               enddo
!
!        ...vertical edge constrained by a rectangular face
            case(25,27,43,46,53)
               nodp = NFACEC(nc)
               call save_min_order(nodp,MAXP*10+nord)
               do k=2,4,2
                  nodp = iabs(NFACE_CONS(k,nc))
                  call save_min_order(nodp,nord)
               enddo
!
!        ...edge constrained by a triangular face
            case(75,76,77)
               nodp = NFACEC(nc)
               call save_min_order(nodp,nord)
               do k=1,3
                  nodp = iabs(NFACE_CONS(k,nc))
                  call save_min_order(nodp,nord)
               enddo
!
!        ...face node constrained by a face
            case(21,22,23,24,31,32,34,35,41,42,44,45,51,52,61,62,71,72,73,74)
               nodp = NFACEC(nc)
               call save_min_order(nodp,nord)
            end select
         endif
      enddo
   enddo
!
!----------------------------------------------------------------------
!  STEP 3: Save the determined order for nodes
!----------------------------------------------------------------------
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
            select case(NODES(nod)%type)
!
!        ...edge node
            case('medg')
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
            case('mdlt')
               call nr_face_sons('mdlt',NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  nods = Son(nod,is)
                  select case(NODES(nods)%type)
                  case('mdlt')
                     nord = min(nord,NODES(nods)%visit)
                  case('mdlq')
                     call decode(NODES(nods)%visit, nordh,nordv)
                     nord = min(nord,nordv)
                  end select
               enddo
!
!           ...modify the node order
               call nodmod(nod, nord)
!
!           ...communicate the new order to the (inactive) sons
               call nr_sons('mdlt',NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  nods = Son(nod,is)
                  select case(NODES(nods)%type)
                  case('medg'); NODES(nods)%order = nord
                  case('mdlt'); NODES(nods)%order = nord
                  case('mdlq'); NODES(nods)%order = nord*10+nord
                  end select
               enddo
!
!        ...rectangular face node
            case('mdlq')
               call decode(nord, nordh,nordv)
               call nr_face_sons('mdlq',NODES(nod)%ref_kind, nrsons)
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
!..reset visitation flags
   call reset_visit
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (RANK .eq. ROOT) write(*,2020) end_time-start_time
 2020 format(' enforce_min: ',f12.5,'  seconds')
!
end subroutine enforce_min_rule
!
!
!
!
subroutine save_min_order(Nod,Nord)
!
   use data_structure3D
   implicit none
!
   integer, intent(in) :: Nod,Nord
!
   integer :: nordh ,nordv
   integer :: nordh1,nordv1
   integer :: nordh2,nordv2
!
   if (NODES(Nod)%visit.eq.0) then
      NODES(nod)%visit= Nord
   else
      select case(NODES(Nod)%type)
      case('medg','mdlt')
         NODES(nod)%visit = min(NODES(nod)%visit,Nord)
      case('mdlq')
         call decode(Nord, nordh1,nordv1)
         call decode(NODES(nod)%visit, nordh2,nordv2)
         nordh = min(nordh1,nordh2); nordv = min(nordv1,nordv2)
         NODES(nod)%visit = nordh*10+nordv
      end select
   endif
!
end subroutine save_min_order
