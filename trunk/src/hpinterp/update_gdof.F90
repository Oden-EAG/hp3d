!-----------------------------------------------------------------------
!
!    routine name       - update_gdof
!
!-----------------------------------------------------------------------
!> @brief Routine updates values of geometry degrees of freedom,
!!        using the GMP parametrizations
!> @date Mar 2023
!-----------------------------------------------------------------------
!
subroutine update_gdof
!
   use constrained_nodes
   use data_structure3D
   use element_data
   use environment, only: QUIET_MODE
   use mpi_wrapper
   use par_mesh   , only: DISTRIBUTED,HOST_MESH
   use GMP
!
   implicit none
!
   integer :: ntype
   character(len=6) :: mdltype
!
!..orientation of element nodes
   integer, dimension(12) :: nedge_orient
   integer, dimension(6)  :: nface_orient
!
!..element nodes
   integer, dimension(27)      :: nodesl, norientl
   integer, dimension(MAXNODM) :: nodm
!
!..order of approximation for an element
   integer, dimension(19) :: norder
!
!..reference coordinates for an element
   real(8), dimension(3,8) :: xsub
!
!..geometry dofs for an element
   real(8), dimension(3,MAXbrickH) :: xnod
!
!..auxiliary variables for timing
   real(8) :: start_time,end_time
!
!..auxiliary variables
   integer :: iel,iv,ie,ifc,ind,iflag,i,k,loc
   integer :: mdle,no,nod,nr_elem_nodes,nrnodm,nr_up_elem
!
!..additional variables for distributed case
   integer :: src,rcv,tag,count,ierr,j_loc,j_glb,j_off,loc_max
   integer :: ndofH,ndofE,ndofV,ndofQ
   logical :: nod_flg
   integer :: nod_cnt(NUM_PROCS)
   integer, allocatable :: nod_loc(:),nod_tmp(:),nod_glb(:),nod_rnk(:)
   real(8), dimension(:,:), pointer :: buf
!
!..use threaded routine; not recommended without proper verification first
   logical, parameter :: USE_THREADED = .false.
   logical, parameter :: opt_blas = .true.
!
!-----------------------------------------------------------------------
!
   if (USE_THREADED) then
      call update_gdof_omp
      return
   endif
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..lower the visitation flag for all nodes
   call reset_visit
!
!-----------------------------------------------------------------------
!
 10 continue
!
!..multiple passes
   do
!
!  ...initiate number of updated elements
      nr_up_elem=0
!
!  ...set flag that indicates if an element was skipped
      nod_flg = .false.
!
!  ...loop through active elements
      do iel=1,NRELES_SUBD
!
         mdle = ELEM_SUBD(iel)
         call find_elem_type(mdle, mdltype)
!
!     ...skip if the element has already been processed
         if (NODES(mdle)%visit.eq.1) cycle
         call refel(mdle, iflag,no,xsub)
!
!     ...determine nodes for the element (active and constrained)
!        and build the data base for the constrained nodes
         call get_connect_info(mdle, nodesl,norientl)
!
!     ...use the info to determine all nodes of the corresponding
!        modified element
         call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
!     ...check if all parent nodes have been updated
         ntype = NODES(mdle)%ntype
         nr_elem_nodes = nvert(ntype)+nedge(ntype)+nface(ntype)+1
         do i=1,nrnodm
            nod = nodm(i)
            call locate(nod,nodesl,nr_elem_nodes, loc)
!
!        ...if not a regular node of the element
            if (loc.eq.0) then
!
!           ...check if the node has been updated
               if (NODES(nod)%visit.eq.0) then
                  nod_flg = .true.
                  cycle
               endif
            endif
         enddo
!
!     ...update the number of processed elements
         nr_up_elem = nr_up_elem+1
!
!-----------------------------------------------------------------------
!
!     ...loop through the element vertex nodes
         do iv=1,nvert(ntype)
            nod = nodesl(iv)
            if (.not.associated(NODES(nod)%dof)) then
               NODES(nod)%visit=1
            endif
            if (NODES(nod)%visit.eq.1) cycle
            call hpvert(iflag,no,xsub(1:3,iv), NODES(nod)%dof%coord)
            NODES(nod)%visit=1
!
!     ...end of loop through element vertices
         enddo
!
!-----------------------------------------------------------------------
!
         call find_orient(mdle, nedge_orient,nface_orient)
         call find_order(mdle, norder)
         call nodcor(mdle, xnod)
!
!     ...loop through element edge nodes
         do ie=1,nedge(ntype)
            ind = nvert(ntype)+ie
            nod = nodesl(ind)
!
!        ...if no gdof, mark as processed
            if (.not.associated(NODES(nod)%dof)) then
               NODES(nod)%visit=1
            elseif (.not.associated(NODES(nod)%dof%coord)) then
               NODES(nod)%visit=1
            endif
            if (NODES(nod)%visit.eq.1) cycle
            if (Is_inactive(nod))            cycle
!
            if (mdltype .ne. 'Linear') then
               call hpedge(mdle,iflag,no,xsub,ntype,             &
                           nedge_orient,nface_orient,norder,ie,  &
                           xnod,NODES(nod)%dof%coord)
            endif
!
            NODES(nod)%visit=1
!
!     ...end of loop through element edges
         enddo
!
!-----------------------------------------------------------------------
!
         call nodcor(mdle, xnod)
!
!     ...loop through element face nodes
         do ifc=1,nface(ntype)
            ind = nvert(ntype)+nedge(ntype)+ifc
            nod = nodesl(ind)
!        ...if no gdof, mark as processed
            if (.not.associated(NODES(nod)%dof)) then
               NODES(nod)%visit=1
            elseif (.not.associated(NODES(nod)%dof%coord)) then
               NODES(nod)%visit=1
            endif
            if (NODES(nod)%visit.eq.1) cycle
            if (Is_inactive(nod))            cycle
!
            if (mdltype .ne. 'Linear') then
               if (opt_blas) then
                  call hpface_opt(mdle,iflag,no,xsub,ntype,             &
                                  nedge_orient,nface_orient,norder,ifc, &
                                  xnod,NODES(nod)%dof%coord)
               else
                  call hpface(mdle,iflag,no,xsub,ntype,             &
                              nedge_orient,nface_orient,norder,ifc, &
                              xnod,NODES(nod)%dof%coord)
               endif
            endif
            NODES(nod)%visit=1
!
!     ...end of loop through element faces
         enddo
!
         NODES(mdle)%visit=1
!
!  ...end of loop through elements
      enddo
!
      if (nr_up_elem.eq.0) exit
!
!..end of infinite loop through multiple passes
   enddo
!
!-----------------------------------------------------------------------
!  Remark: the following node data exchange is needed to take care of
!          constrained nodes at the subdomain interface where the parent
!          node of the other subdomain must be updated first, before the
!          constrained nodes can be updated.
!
!                        subdomain interface
!                                 .
!                                 .
!                     *-----------*-----------*
!                     |           |           |
!                     |           |           |
!                     *-----------*           |
!                     |           |           |
!                     |           |           |
!                     *-----------*-----------*
!                                 .
!                                 .
!-----------------------------------------------------------------------
!
!..skip node data exchange if this is not a distributed run
   if ((.not. DISTRIBUTED) .or. HOST_MESH) then
!  ...consistency check
      if (nod_flg) then
         write(*,*) 'update_gdof: mesh is not distributed (or host mesh), but flag = ',nod_flg
         stop
      endif
      goto 70
   endif
!
   !write(*,4010) RANK,nod_flg
   4010 format('update_gdof: [',I4,'], nod_flag=',L2)
!
!..skip creating node list if no DOFs are needed from other subdomains
   j_loc = 0
   if (.not. nod_flg) goto 60
!
!..fill local node list
   loc_max = 1000; allocate(nod_loc(loc_max))
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!  ...skip if the element has already been processed
      if (NODES(mdle)%visit.eq.1) cycle
!  ...determine nodes for the modified element
      call get_connect_info(mdle, nodesl,norientl)
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
      nr_elem_nodes = nvert(ntype)+nedge(ntype)+nface(ntype)+1
!  ...iterate through the nodes of the modified element
      do i=1,nrnodm
         nod = nodm(i)
         call locate(nod,nodesl,nr_elem_nodes, loc)
!     ...check if the node has been updated
         if (loc.eq.0 .and. NODES(nod)%visit.eq.0) then
            NODES(nod)%visit = -1
            if(j_loc .ge. loc_max) then
               loc_max = loc_max*2
               allocate(nod_tmp(loc_max))
               nod_tmp(1:j_loc) = nod_loc(1:j_loc)
               call move_alloc(nod_tmp, nod_loc)
            endif
            j_loc = j_loc+1
            nod_loc(j_loc) = nod
         endif
      enddo
   enddo
!
!..consistency check
   if (j_loc .eq. 0) then
      write(*,*) 'update_gdof: j_loc = 0, but flag = ',nod_flg
      stop
   endif
!
   !write(*,5010) RANK,j_loc
 5010 format('update_gdof: [',I4,'], jloc=',I6)
!
 60 continue
!
!..collect local counts from every processor
   nod_cnt(1:NUM_PROCS) = 0; nod_cnt(RANK+1) = j_loc
   count = NUM_PROCS
   call MPI_ALLREDUCE(MPI_IN_PLACE,nod_cnt,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!..compute global count and offset
   j_off = 0
   if (RANK > 0) j_off = sum(nod_cnt(1:RANK))
   j_glb = sum(nod_cnt(1:NUM_PROCS))
!
!..exit if nobody needs DOFs from other subdomains
   if (j_glb .eq. 0) goto 70
!
!..fill global node list with local nodes
   allocate(nod_glb(j_glb)); nod_glb = 0
   if (nod_flg) then
      nod_glb(j_off+1:j_off+j_loc) = nod_loc(1:j_loc)
      deallocate(nod_loc)
   endif
!
!..collect local node list from every processor
   count = j_glb
   call MPI_ALLREDUCE(MPI_IN_PLACE,nod_glb,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!..iterate through the node list and check whether one can supply node data
   allocate(nod_rnk(j_glb)); nod_rnk = NUM_PROCS
   do i=1,j_glb
      nod = nod_glb(i)
      if (NODES(nod)%visit.eq.1) nod_rnk(i) = RANK
   enddo
!
!..collect node supplier list
   count = j_glb
   call MPI_ALLREDUCE(MPI_IN_PLACE,nod_rnk,count,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
!
!..iterate through node supplier list, and exchange node data
   rcv = 0; k = nod_cnt(rcv+1)
   do i=1,j_glb
!  ...check whether there is no supplier
      if (nod_rnk(i) .eq. NUM_PROCS) cycle
!  ...determine supplier, receiver, and dof count
      src = nod_rnk(i)
      tag = nod_glb(i)
      do while (i > k)
         rcv = rcv+1
         k = k+nod_cnt(rcv+1)
      enddo
      nod = nod_glb(i)
      call find_ndof(nod, ndofH,ndofE,ndofV,ndofQ)
      if (ndofH .eq. 0) then
         NODES(nod)%visit = 1
         cycle
      endif
      count = 3*ndofH
!  ...check whether supplier
      if (src .eq. RANK) then
         if (.not. associated(NODES(nod)%dof)) then
            write(*,6010) RANK,'SRC: DOF NOT ASSOCIATED  ',nod,S_Type(NODES(nod)%ntype)
            stop
         endif
         if (.not. associated(NODES(nod)%dof%coord)) then
            write(*,6010) RANK,'SRC: COORD NOT ASSOCIATED  ',nod,S_Type(NODES(nod)%ntype)
            stop
         endif
         !write(*,6010) RANK,'SENDING  ',nod,S_Type(NODES(nod)%ntype)
         buf => NODES(nod)%dof%coord
         call MPI_SEND(buf,count,MPI_REAL8,rcv,tag,MPI_COMM_WORLD,ierr)
!  ...check whether receiver
      elseif (rcv .eq. RANK) then
         if (.not. associated(NODES(nod)%dof)) then
            write(*,6010) RANK,'RCV: DOF NOT ASSOCIATED  ',nod,S_Type(NODES(nod)%ntype)
            stop
         endif
         if (.not. associated(NODES(nod)%dof%coord)) then
            write(*,6010) RANK,'RCV: COORD NOT ASSOCIATED  ',nod,S_Type(NODES(nod)%ntype)
            stop
         endif
         !write(*,6010) RANK,'RECEIVING',nod,S_Type(NODES(nod)%ntype)
         buf => NODES(nod)%dof%coord
         call MPI_RECV(buf,count,MPI_REAL8,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         NODES(nod)%visit = 1
      endif
 6010 format('update_gdof: [',I4,'] ',A,' nod=',I6,', type=',A)
   enddo
!
   if (nod_flg) then
      do i=1,j_loc
         nod = nod_glb(j_off+i)
!     ...reset visitation flag if node was not supplied
         if (NODES(nod)%visit.eq.-1) then
            NODES(nod)%visit = 0
         endif
      enddo
   endif
!
   deallocate(nod_rnk,nod_glb)
!
!..wait for non-blocking send/recv to finish (with MPI_ISEND/MPI_IRECV)
   !call MPI_WAITALL()
!
!..go back to first loop and process elements if necessary
   if (nod_flg) then
      goto 10
   else
      goto 60
   endif
!
 70 continue
!
   if (allocated(nod_loc)) deallocate(nod_loc)
   if (allocated(nod_glb)) deallocate(nod_glb)
   if (allocated(nod_rnk)) deallocate(nod_rnk)
!
!-----------------------------------------------------------------------
!
!..update the middle node
!
!$OMP PARALLEL DO                                           &
!$OMP PRIVATE(iflag,mdle,mdltype,nedge_orient,nface_orient, &
!$OMP         no,norder,ntype,xnod,xsub)                    &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!  ...cycle if middle node does not have H1 DOFs (low order)
      if (.not.associated(NODES(mdle)%dof%coord)) cycle
      call find_elem_type(mdle, mdltype)
      if (mdltype .ne. 'Linear') then
         ntype = NODES(mdle)%ntype
         call refel(mdle, iflag,no,xsub)
         call find_orient(mdle, nedge_orient,nface_orient)
         call find_order(mdle, norder)
         call nodcor(mdle, xnod)
         if (opt_blas) then
            call hpmdle_opt(mdle,iflag,no,xsub,ntype,          &
                            nedge_orient,nface_orient,norder,  &
                            xnod,NODES(mdle)%dof%coord)
         else
            call hpmdle(mdle,iflag,no,xsub,ntype,          &
                        nedge_orient,nface_orient,norder,  &
                        xnod,NODES(mdle)%dof%coord)
         endif
      endif
   enddo
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,8004) end_time-start_time
 8004 format(' update_gdof: ',f12.5,'  seconds')
   endif
!
end subroutine update_gdof
!
!
!
!-----------------------------------------------------------------------
!
!    routine name       - update_gdof_omp
!
!-----------------------------------------------------------------------
!> @brief OpenMP parallel version of update_gdof
!> @date Mar 2023
!-----------------------------------------------------------------------
!
subroutine update_gdof_omp
!
   use constrained_nodes
   use data_structure3D
   use element_data
   use environment, only: QUIET_MODE
   use mpi_wrapper
   use par_mesh   , only: DISTRIBUTED,HOST_MESH
   use GMP
!
   implicit none
!
   integer :: ntype
   character(len=6) :: mdltype
!
!..orientation of element nodes
   integer, dimension(12) :: nedge_orient
   integer, dimension(6)  :: nface_orient
!
!..element nodes
   integer, dimension(27)      :: nodesl, norientl
   integer, dimension(MAXNODM) :: nodm
!
!..order of approximation for an element
   integer, dimension(19) :: norder
!
!..reference coordinates for an element
   real(8), dimension(3,8) :: xsub
!
!..geometry dofs for an element
   real(8), dimension(3,MAXbrickH) :: xnod
!
!..temporary for OpenMP parallelism
   real(8) :: coord(3,MAXquadH)
!
!..auxiliary variables for timing
   real(8) :: start_time,end_time
!
!..auxiliary variables
   integer :: iel,iv,ie,ifc,ind,iflag,i,k,loc
   integer :: mdle,no,nod,nr_elem_nodes,nrnodm,nr_up_elem
!
!..additional variables for distributed case
   integer :: src,rcv,tag,count,ierr,j_loc,j_glb,j_off,loc_max
   integer :: ndofH,ndofE,ndofV,ndofQ
   logical :: nod_flg
   integer :: nod_cnt(NUM_PROCS)
   integer, allocatable :: nod_loc(:),nod_tmp(:),nod_glb(:),nod_rnk(:)
   real(8), dimension(:,:), pointer :: buf
!
!..OpenMP parallelization
   integer, allocatable :: nodes_elem(:,:), nodes_irreg(:,:), nrnodes_irreg(:)
   integer :: nodesi(27)
   integer :: j, sz, nrnodi
!
!-----------------------------------------------------------------------
!
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   start_time = MPI_Wtime()
!
   call reset_visit
!
!-----------------------------------------------------------------------
!
   allocate(nodes_elem(27,NRELES_SUBD))
   allocate(nodes_irreg(27,NRELES_SUBD))
   allocate(nrnodes_irreg(NRELES_SUBD))
!
!-----------------------------------------------------------------------
!
! Get nodes for each element in subd
!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(iel,mdle,nodesl,norientl,nodm,nrnodm,   &
!$OMP            nr_elem_nodes,loc,nod,i,j,ntype)        &
!$OMP    SCHEDULE(STATIC) REDUCTION(+:nr_up_elem)
!..loop through active elements
   do iel=1,NRELES_SUBD
!
      mdle = ELEM_SUBD(iel)
!
!  ...get unconstrained nodes
      call get_connect_info(mdle, nodesl,norientl)
!  ...get constrained nodes
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
!  ...check if constrained nodes are regular (negative for irregular)
      ntype = NODES(mdle)%ntype
      nr_elem_nodes = nvert(ntype)+nedge(ntype)+nface(ntype)+1
      do i=1,nr_elem_nodes
         nod = nodesl(i)
         call locate(nod,nodm,nrnodm, loc)
!
!     ...if not a regular node of the element, negate
         if (loc.eq.0) then
            nodes_elem(i,iel) = -nod
         else
            nodes_elem(i,iel) = nod
            nodm(loc) = 0
         endif
      enddo
!
      j = 0
      do i=1,nrnodm
         nod = nodm(i)
         if (nod.ne.0) then
            j = j + 1
            nodes_irreg(j,iel) = nod
         endif
      enddo
      nrnodes_irreg(iel) = j
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
 10 continue
!
!..multiple passes
   do
!
!  ...initiate number of updated elements
      nr_up_elem=0
!
!  ...set flag that indicates if an element was skipped
      nod_flg = .false.
!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(iel,mdle,iflag,no,xsub,nrnodi,nodesi,nodesl,sz, &
!$OMP            nr_elem_nodes,nface_orient,iv,ie,ifc,nod,coord, &
!$OMP            ntype,norder,ind,mdltype,nedge_orient,xnod)     &
!$OMP    SCHEDULE(DYNAMIC) REDUCTION(+:nr_up_elem)
!  ...loop through active elements
      do iel=1,NRELES_SUBD
!
         mdle = ELEM_SUBD(iel)
         call find_elem_type(mdle, mdltype)
!
!     ...skip if the element has already been processed
         if (NODES(mdle)%visit.eq.1) cycle
         call refel(mdle, iflag,no,xsub)
!
!     ...check if all parent nodes have been updated
         ntype = NODES(mdle)%ntype
         nr_elem_nodes = nvert(ntype)+nedge(ntype)+nface(ntype)+1
!
         nrnodi = nrnodes_irreg(iel)
         nodesi(1:nrnodi) = nodes_irreg(1:nrnodi,iel)
         do i=1,nrnodi
            nod = nodesi(i)
!
!        ...check if the node has been updated
            if (NODES(nod)%visit.eq.0) then
               nod_flg = .true.
               cycle
            endif
         enddo
!
         nodesl(1:27) = nodes_elem(1:27,iel)
!
!     ...update the number of processed elements
         nr_up_elem = nr_up_elem+1
!
!-----------------------------------------------------------------------
!
!     ...loop through the element vertex nodes
         do iv=1,nvert(ntype)
            nod = abs(nodesl(iv))
            if (.not.associated(NODES(nod)%dof)) then
               NODES(nod)%visit=1
            endif
            if (NODES(nod)%visit.eq.1) cycle
!
            sz = size(NODES(nod)%dof%coord,2)
!
            call hpvert(iflag,no,xsub(1:3,iv), coord(1:3,1:sz))
!$OMP CRITICAL
            if (NODES(nod)%visit.eq.0) then
               NODES(nod)%dof%coord = coord(1:3,1:sz)
               NODES(nod)%visit=1
            endif
!$OMP END CRITICAL
!
!     ...end of loop through element vertices
         enddo
!
!-----------------------------------------------------------------------
!
         call find_orient(mdle, nedge_orient,nface_orient)
         call find_order(mdle, norder)
         call nodcor(mdle, xnod)
!
!     ...loop through element edge nodes
         do ie=1,nedge(ntype)
            ind = nvert(ntype)+ie
            nod = iabs(nodesl(ind))
!
!        ...if no gdof, mark as processed
            if (.not.associated(NODES(nod)%dof)) then
               NODES(nod)%visit=1
            elseif (.not.associated(NODES(nod)%dof%coord)) then
               NODES(nod)%visit=1
            endif
            if (NODES(nod)%visit.eq.1) cycle
            if (Is_inactive(nod))            cycle
!
            sz = size(NODES(nod)%dof%coord,2)
!
            if (mdltype .ne. 'Linear') then
               call hpedge(mdle,iflag,no,xsub,ntype,             &
                           nedge_orient,nface_orient,norder,ie,  &
                           xnod,coord(1:3,1:sz))
!
!$OMP CRITICAL
               if (NODES(nod)%visit.eq.0) then
                  NODES(nod)%dof%coord = coord(1:3,1:sz)
                  NODES(nod)%visit=1
               endif
!$OMP END CRITICAL
!
            else
               NODES(nod)%visit=1
            endif
!
!     ...end of loop through element edges
         enddo
!
!-----------------------------------------------------------------------
!
         call nodcor(mdle, xnod)
!
!     ...loop through element face nodes
         do ifc=1,nface(ntype)
            ind = nvert(ntype)+nedge(ntype)+ifc
            nod = iabs(nodesl(ind))
!        ...if no gdof, mark as processed
            if (.not.associated(NODES(nod)%dof)) then
               NODES(nod)%visit=1
            elseif (.not.associated(NODES(nod)%dof%coord)) then
               NODES(nod)%visit=1
            endif
            if (NODES(nod)%visit.eq.1) cycle
            if (Is_inactive(nod))            cycle
!
            sz = size(NODES(nod)%dof%coord,2)
!
            if (mdltype .ne. 'Linear') then
               call hpface_opt(mdle,iflag,no,xsub,ntype,             &
                               nedge_orient,nface_orient,norder,ifc, &
                               xnod,coord(1:3,1:sz))
!
!$OMP CRITICAL
               if (NODES(nod)%visit.eq.0) then
                  NODES(nod)%dof%coord = coord(1:3,1:sz)
                  NODES(nod)%visit=1
               endif
!$OMP END CRITICAL
!
            else
               NODES(nod)%visit=1
            endif
!
!     ...end of loop through element faces
         enddo
!
         NODES(mdle)%visit=1
!
!  ...end of loop through elements
      enddo
!$OMP END DO
!$OMP END PARALLEL
!
      if (nr_up_elem.eq.0) exit
!
!..end of infinite loop through multiple passes
   enddo
!
!-----------------------------------------------------------------------
!  Remark: the following node data exchange is needed to take care of
!          constrained nodes at the subdomain interface where the parent
!          node of the other subdomain must be updated first, before the
!          constrained nodes can be updated.
!
!                        subdomain interface
!                                 .
!                                 .
!                     *-----------*-----------*
!                     |           |           |
!                     |           |           |
!                     *-----------*           |
!                     |           |           |
!                     |           |           |
!                     *-----------*-----------*
!                                 .
!                                 .
!-----------------------------------------------------------------------
!
!..skip node data exchange if this is not a distributed run
   if ((.not. DISTRIBUTED) .or. HOST_MESH) then
!  ...consistency check
      if (nod_flg) then
         write(*,*) 'update_gdof: mesh is not distributed (or host mesh), but flag = ',nod_flg
         stop
      endif
      goto 70
   endif
!
   !write(*,4010) RANK,nod_flg
   4010 format('update_gdof: [',I4,'], nod_flag=',L2)
!
!..skip creating node list if no DOFs are needed from other subdomains
   j_loc = 0
   if (.not. nod_flg) goto 60
!
!..fill local node list
   loc_max = 1000; allocate(nod_loc(loc_max))
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!  ...skip if the element has already been processed
      if (NODES(mdle)%visit.eq.1) cycle
!
      ntype = NODES(mdle)%ntype
      nr_elem_nodes = nvert(ntype)+nedge(ntype)+nface(ntype)+1
!
!  ...check if irregular nodes have been updated; add to list if not
      nrnodi = nrnodes_irreg(iel)
      nodesi(1:nrnodi) = nodes_irreg(1:nrnodi,iel)
      do i=1,nrnodi
         nod = nodesi(i)
!
!     ...check if the node has been updated
         if (NODES(nod)%visit.eq.0) then
            NODES(nod)%visit = -1
            if(j_loc .ge. loc_max) then
               loc_max = loc_max*2
               allocate(nod_tmp(loc_max))
               nod_tmp(1:j_loc) = nod_loc(1:j_loc)
               call move_alloc(nod_tmp, nod_loc)
            endif
            j_loc = j_loc+1
            nod_loc(j_loc) = nod
         endif
      enddo
   enddo
!
!..consistency check
   if (j_loc .eq. 0) then
      write(*,*) 'update_gdof: j_loc = 0, but flag = ',nod_flg
      stop
   endif
!
!  write(*,5010) RANK,j_loc
 5010 format('update_gdof: [',I4,'], jloc=',I6)
!
 60 continue
!
!..collect local counts from every processor
   nod_cnt(1:NUM_PROCS) = 0; nod_cnt(RANK+1) = j_loc
   count = NUM_PROCS
   call MPI_ALLREDUCE(MPI_IN_PLACE,nod_cnt,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!..compute global count and offset
   j_off = 0
   if (RANK > 0) j_off = sum(nod_cnt(1:RANK))
   j_glb = sum(nod_cnt(1:NUM_PROCS))
!
!..exit if nobody needs DOFs from other subdomains
   if (j_glb .eq. 0) goto 70
!
!..fill global node list with local nodes
   allocate(nod_glb(j_glb)); nod_glb = 0
   if (nod_flg) then
      nod_glb(j_off+1:j_off+j_loc) = nod_loc(1:j_loc)
      deallocate(nod_loc)
   endif
!
!..collect local node list from every processor
   count = j_glb
   call MPI_ALLREDUCE(MPI_IN_PLACE,nod_glb,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!..iterate through the node list and check whether one can supply node data
   allocate(nod_rnk(j_glb)); nod_rnk = NUM_PROCS
   do i=1,j_glb
      nod = nod_glb(i)
      if (NODES(nod)%visit .eq. 1) nod_rnk(i) = RANK
   enddo
!
!..collect node supplier list
   count = j_glb
   call MPI_ALLREDUCE(MPI_IN_PLACE,nod_rnk,count,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
!
!..iterate through node supplier list, and exchange node data
   rcv = 0; k = nod_cnt(rcv+1)
   do i=1,j_glb
!  ...check whether there is no supplier
      if (nod_rnk(i) .eq. NUM_PROCS) cycle
!  ...determine supplier, receiver, and dof count
      src = nod_rnk(i)
      tag = mod(nod_glb(i),2000000)
      do while (i > k)
         rcv = rcv+1
         k = k+nod_cnt(rcv+1)
      enddo
      nod = nod_glb(i)
      call find_ndof(nod, ndofH,ndofE,ndofV,ndofQ)
      if (ndofH .eq. 0) then
         NODES(nod)%visit = 1
         cycle
      endif
      count = 3*ndofH
!  ...check whether supplier
      if (src .eq. RANK) then
         if (.not. associated(NODES(nod)%dof)) then
            write(*,6010) RANK,'SRC: DOF NOT ASSOCIATED  ',nod,S_Type(NODES(nod)%ntype)
            stop
         endif
         if (.not. associated(NODES(nod)%dof%coord)) then
            write(*,6010) RANK,'SRC: COORD NOT ASSOCIATED  ',nod,S_Type(NODES(nod)%ntype)
            stop
         endif
         !write(*,6010) RANK,'SENDING  ',nod,S_Type(NODES(nod)%ntype)
         buf => NODES(nod)%dof%coord
         call MPI_SEND(buf,count,MPI_REAL8,rcv,tag,MPI_COMM_WORLD,ierr)
!  ...check whether receiver
      elseif (rcv .eq. RANK) then
         if (.not. associated(NODES(nod)%dof)) then
            write(*,6010) RANK,'RCV: DOF NOT ASSOCIATED  ',nod,S_Type(NODES(nod)%ntype)
            stop
         endif
         if (.not. associated(NODES(nod)%dof%coord)) then
            write(*,6010) RANK,'RCV: COORD NOT ASSOCIATED  ',nod,S_Type(NODES(nod)%ntype)
            stop
         endif
         !write(*,6010) RANK,'RECEIVING',nod,S_Type(NODES(nod)%ntype)
         buf => NODES(nod)%dof%coord
         call MPI_RECV(buf,count,MPI_REAL8,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         NODES(nod)%visit = 1
      endif
 6010 format('update_gdof: [',I4,'] ',A,' nod=',I6,', type=',A)
   enddo
!
   if (nod_flg) then
      do i=1,j_loc
         nod = nod_glb(j_off+i)
!     ...reset visitation flag if node was not supplied
         if (NODES(nod)%visit.eq.-1) then
            NODES(nod)%visit = 0
         endif
      enddo
   endif
!
   deallocate(nod_rnk,nod_glb)
!
!..wait for non-blocking send/recv to finish (with MPI_ISEND/MPI_IRECV)
   !call MPI_WAITALL()
!
!..go back to first loop and process elements if necessary
   if (nod_flg) then
      goto 10
   else
      goto 60
   endif
!
 70 continue
!
   if (allocated(nod_loc)) deallocate(nod_loc)
   if (allocated(nod_glb)) deallocate(nod_glb)
   if (allocated(nod_rnk)) deallocate(nod_rnk)
   deallocate(nodes_elem)
   deallocate(nodes_irreg)
   deallocate(nrnodes_irreg)
!
!-----------------------------------------------------------------------
!
!..update the middle node
!
!$OMP PARALLEL DO                                           &
!$OMP PRIVATE(iflag,mdle,mdltype,nedge_orient,nface_orient, &
!$OMP         no,norder,ntype,xnod,xsub)                    &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!  ...cycle if middle node does not have H1 DOFs (low order)
      if (.not.associated(NODES(mdle)%dof%coord)) cycle
      call find_elem_type(mdle, mdltype)
!
      if (mdltype .ne. 'Linear') then
         ntype = NODES(mdle)%ntype
         call refel(mdle, iflag,no,xsub)
         call find_orient(mdle, nedge_orient,nface_orient)
         call find_order(mdle, norder)
         call nodcor(mdle, xnod)
         call hpmdle_opt(mdle,iflag,no,xsub,ntype,          &
                         nedge_orient,nface_orient,norder,  &
                         xnod,NODES(mdle)%dof%coord)
      endif
   enddo
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      end_time = MPI_Wtime()
      write(*,8004) end_time-start_time
 8004 format(' update_gdof: ',f12.5,'  seconds')
   endif
!
end subroutine update_gdof_omp

