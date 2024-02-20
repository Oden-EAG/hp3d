#include "typedefs.h"
!-----------------------------------------------------------------------
!
!    routine name       - update_Ddof
!
!-----------------------------------------------------------------------
!> @brief Routine updates values of solution degrees of freedom
!!        for Dirichlet nodes
!> @date Sep 2023
!-----------------------------------------------------------------------
!
subroutine update_Ddof
!
   use data_structure3D
   use environment, only: QUIET_MODE
   use par_mesh   , only: DISTRIBUTED
   use MPI        , only: MPI_COMM_WORLD,MPI_INTEGER,MPI_REAL8,MPI_COMPLEX16, &
                          MPI_SUM,MPI_MIN,MPI_IN_PLACE,MPI_STATUS_IGNORE,     &
                          MPI_BCAST,MPI_REDUCE,MPI_ALLREDUCE,MPI_Wtime
   use mpi_param  , only: RANK,ROOT,NUM_PROCS
   use par_mesh   , only: DISTRIBUTED,HOST_MESH
!
   implicit none
!
   integer :: ntype
!
!..orientation of element nodes
   integer, dimension(12)      :: nedge_orient
   integer, dimension(6)       :: nface_orient
!
!..element nodes
   integer, dimension(27)      :: nodesl, norientl
   integer, dimension(MAXNODM) :: nodm
!
!..order of approximation for an element
   integer, dimension(19)      :: norder
!
!..reference coordinates for an element
   real(8), dimension(3,8)     :: xsub
!
!..solution dof for an element
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!..auxiliary variables for timing
   real(8) :: start_time,end_time
!
!..auxiliary variables
   integer :: iel,iv,ie,ifc,ind,iflag
   integer :: mdle,no,nod
   integer :: i,k,loc,nr_elem_nodes,nrnodm,nr_up_elem
!
!..additional variables for distributed case
   integer :: src,rcv,tag,count,ierr,j_loc,j_glb,j_off,loc_max
   integer :: ndofH,ndofE,ndofV,ndofQ,icase,nvarH,nvarE,nvarV
   logical :: nod_flg
   integer :: nod_cnt(NUM_PROCS)
   integer, allocatable :: nod_loc(:),nod_tmp(:),nod_glb(:),nod_rnk(:)
   VTYPE  , dimension(:,:), pointer :: buf
!
!..WARNING: "dirichlet" and other user-supplied routines must be thread-
!           safe to use this routine; not recommended without proper
!           verification
   logical :: USE_THREADED = .false.
!
#if DEBUG_MODE
   integer :: iprint
   iprint=0
#endif
!
!-----------------------------------------------------------------------
!
   if (USE_THREADED) then
      call update_Ddof_omp
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
      do 100 iel=1,NRELES_SUBD
!
         mdle = ELEM_SUBD(iel)
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
               if (NODES(nod)%visit.eq.0 .and. is_Dirichlet(nod)) then
                  nod_flg = .true.
                  goto 100
               endif
            endif
         enddo
!
!     ...update the number of processed elements
         nr_up_elem = nr_up_elem+1
!
!-----------------------------------------------------------------------
!     ...Step 1: Update   V E R T   dof for Dirichlet nodes
         do iv=1,nvert(ntype)
            nod = nodesl(iv)
            if (.not.associated(NODES(nod)%dof))       cycle
            if (.not.associated(NODES(nod)%dof%zdofH)) cycle
            if (NODES(nod)%visit.eq.1)                 cycle
            if (is_Dirichlet_attr(nod,CONTIN)) then
#if DEBUG_MODE
               if (iprint.eq.1) write(*,7010) mdle,iv,nod
          7010 format('update_Ddof: CALLING dhpvert FOR mdle,iv,nod = ',i8,i2,i8)
#endif
               call dhpvert(mdle,iflag,no,xsub(1:3,iv),NODES(nod)%case, &
                            NODES(nod)%bcond, NODES(nod)%dof%zdofH(:,:,N_COMS))
               NODES(nod)%visit=1
            endif
         enddo
!
!-----------------------------------------------------------------------
!
!     ...Step 2: Update   E D G E   dof for Dirichlet nodes
         call find_orient(mdle, nedge_orient,nface_orient)
         call find_order(mdle, norder)
!
!     ...compute solution dofs (need for H1 update)
         call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!
         do ie=1,nedge(ntype)
            ind = nvert(ntype)+ie
            nod = nodesl(ind)
            if (.not.associated(NODES(nod)%dof)) cycle
            if (NODES(nod)%visit.eq.1)           cycle
            if (is_Dirichlet_attr(nod,CONTIN)) then
!           ...update H1 Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofH)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7020) mdle,ie,nod
             7020 format('update_Ddof: CALLING dhpedgeH FOR mdle,ie,nod = ',i8,i2,i8)
#endif
                  call dhpedgeH(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ie,    &
                                zdofH, NODES(nod)%dof%zdofH(:,:,N_COMS))
               endif
               NODES(nod)%visit=1
            endif
            if (is_Dirichlet_attr(nod,TANGEN)) then
!           ...update H(curl) Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofE)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7030) mdle,ie,nod
             7030 format('update_Ddof: CALLING dhpedgeE FOR mdle,ie,nod = ',i8,i2,i8)
#endif
                  call dhpedgeE(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ie,    &
                                NODES(nod)%dof%zdofE(:,:,N_COMS))
               endif
               NODES(nod)%visit=1
            endif
         enddo
!
!-----------------------------------------------------------------------
!  ...Step 3: Update   F A C E   dof for Dirichlet nodes
!
!     ...compute solution dofs (needed for H1 and Hcurl update)
         call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!
!     ...loop over faces
         do ifc=1,nface(ntype)
!        ...get local node number
            ind = nvert(ntype)+nedge(ntype)+ifc
!        ...get global node number
            nod = nodesl(ind)
            if (.not.associated(NODES(nod)%dof)) cycle
            if (NODES(nod)%visit.eq.1)           cycle
            if (is_Dirichlet_attr(nod,CONTIN)) then
!           ...update H1 Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofH)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7040) mdle,ifc,nod
             7040 format('update_Ddof: CALLING dhpfaceH FOR mdle,ifc,nod = ',i8,i2,i8)
#endif
                  call dhpfaceH(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ifc,   &
                                zdofH, NODES(nod)%dof%zdofH(:,:,N_COMS))
               endif
               NODES(nod)%visit=1
            endif
            if (is_Dirichlet_attr(nod,TANGEN)) then
!           ...update H(curl) Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofE)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7050) mdle,ifc,nod
             7050 format('update_Ddof: CALLING dhpfaceE FOR mdle,ifc,nod = ',i8,i2,i8)
#endif
                  call dhpfaceE(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ifc,   &
                                zdofE, NODES(nod)%dof%zdofE(:,:,N_COMS))
               endif
               NODES(nod)%visit=1
            endif
            if (is_Dirichlet_attr(nod,NORMAL)) then
!           ...update H(div) Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofV)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7060) mdle,ifc,nod
             7060 format('update_Ddof: CALLING dhpfaceV FOR mdle,ifc,nod = ',i8,i2,i8)
#endif
                  call dhpfaceV(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ifc,   &
                                NODES(nod)%dof%zdofV(:,:,N_COMS))
               endif
               NODES(nod)%visit=1
            endif
!     ...end of loop over faces
         enddo
!
         NODES(mdle)%visit=1
!
!  ...end of loop through elements
 100  continue
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
         write(*,*) 'update_Ddof: mesh is not distributed (or host mesh), but flag = ',nod_flg
         stop
      endif
      goto 70
   endif
!
   !write(*,4010) RANK,nod_flg
 4010 format('update_Ddof: [',I4,'], nod_flag=',L2)
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
         if (loc.eq.0 .and. NODES(nod)%visit.eq.0  &
                      .and. is_Dirichlet(nod) ) then
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
      write(*,*) 'update_Ddof: j_loc = 0, but flag = ',nod_flg
      stop
   endif
!
   !write(*,5010) RANK,j_loc
   5010 format('update_Ddof: [',I4,'], jloc=',I6)
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
      icase = NODES(Nod)%case
      nvarH = NREQNH(icase)*NRRHS
      nvarE = NREQNE(icase)*NRRHS
      nvarV = NREQNV(icase)*NRRHS
!  ...check whether supplier
      if (src .eq. RANK) then
         !write(*,6010) RANK,'SENDING  ',nod,S_Type(NODES(nod)%ntype)
         if (nvarH > 0 .and. ndofH > 0) then
            count = ndofH*nvarH; buf => NODES(nod)%dof%zdofH(:,:,N_COMS)
            call MPI_SEND(buf,count,MPI_VTYPE,rcv,tag,MPI_COMM_WORLD,ierr)
         endif
         if (nvarE > 0 .and. ndofE > 0) then
            count = ndofE*nvarE; buf => NODES(nod)%dof%zdofE(:,:,N_COMS)
            call MPI_SEND(buf,count,MPI_VTYPE,rcv,tag,MPI_COMM_WORLD,ierr)
         endif
         if (nvarV > 0 .and. ndofV > 0) then
            count = ndofV*nvarV; buf => NODES(nod)%dof%zdofV(:,:,N_COMS)
            call MPI_SEND(buf,count,MPI_VTYPE,rcv,tag,MPI_COMM_WORLD,ierr)
         endif
!  ...check whether receiver
      elseif (rcv .eq. RANK) then
         !write(*,6010) RANK,'RECEIVING',nod,S_Type(NODES(nod)%ntype)
         if (nvarH > 0 .and. ndofH > 0) then
            count = ndofH*nvarH; buf => NODES(nod)%dof%zdofH(:,:,N_COMS)
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         endif
         if (nvarE > 0 .and. ndofE > 0) then
            count = ndofE*nvarE; buf => NODES(nod)%dof%zdofE(:,:,N_COMS)
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         endif
         if (nvarV > 0 .and. ndofV > 0) then
            count = ndofV*nvarV; buf => NODES(nod)%dof%zdofV(:,:,N_COMS)
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         endif
         NODES(nod)%visit = 1
      endif
 6010 format('update_Ddof: [',I4,'] ',A,' nod=',I6,', type=',A)
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
!..wait for non-blocking send/recv to finish (use MPI_ISEND/MPI_IRECV)
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
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,8010) end_time-start_time
 8010 format(' update_Ddof: ',f12.5,'  seconds',/)
   endif
!
end subroutine update_Ddof
!
!
!
!-----------------------------------------------------------------------
!
!    routine name       - update_Ddof_omp
!
!-----------------------------------------------------------------------
!> @brief OpenMP parallel version of update_Ddof
!> @warning  "dirichlet" and other user-supplied routines must be
!!           threadsafe to use this version; not recommended without
!!           proper verification.
!> @date Sep 2023
!-----------------------------------------------------------------------
!
subroutine update_Ddof_omp
!
   use data_structure3D
   use environment, only: QUIET_MODE
   use par_mesh   , only: DISTRIBUTED
   use MPI        , only: MPI_COMM_WORLD,MPI_INTEGER,MPI_REAL8,MPI_COMPLEX16, &
                          MPI_SUM,MPI_MIN,MPI_IN_PLACE,MPI_STATUS_IGNORE
   use mpi_param  , only: RANK,ROOT,NUM_PROCS
   use par_mesh   , only: DISTRIBUTED,HOST_MESH
!
   implicit none
!
   integer :: ntype
!
!..orientation of element nodes
   integer, dimension(12)      :: nedge_orient
   integer, dimension(6)       :: nface_orient
!
!..element nodes
   integer, dimension(27)      :: nodesl, norientl
   integer, dimension(MAXNODM) :: nodm
!
!..order of approximation for an element
   integer, dimension(19)      :: norder
!
!..reference coordinates for an element
   real(8), dimension(3,8)     :: xsub
!
!..solution dof for an element
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!..auxiliary variables for timing
   real(8) :: MPI_Wtime,start_time,end_time
!
!..auxiliary variables
   integer :: iel,iv,ie,ifc,ind,iflag
   integer :: mdle,no,nod
   integer :: i,k,loc,nr_elem_nodes,nrnodm,nr_up_elem
!
!..additional variables for distributed case
   integer :: src,rcv,tag,count,ierr,j_loc,j_glb,j_off,loc_max
   integer :: ndofH,ndofE,ndofV,ndofQ,icase,nvarH,nvarE,nvarV
   logical :: nod_flg
   integer :: nod_cnt(NUM_PROCS)
   integer, allocatable :: nod_loc(:),nod_tmp(:),nod_glb(:),nod_rnk(:)
   VTYPE  , dimension(:,:), pointer :: buf
!
!..OpenMP parallelization
   integer, allocatable :: nodes_elem(:,:), nodes_irreg(:,:), nrnodes_irreg(:)
   integer :: nodesi(27)
   integer :: j, nvar, ndof, nrnodi
!
   VTYPE :: tempH(MAXEQNH,MAXquadH)
   VTYPE :: tempE(MAXEQNE,MAXquadE)
   VTYPE :: tempV(MAXEQNV,MAXquadV)
!
#if DEBUG_MODE
   integer :: iprint = 0
#endif
!
!-----------------------------------------------------------------------
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..lower the visitation flag for all nodes
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
!  ...loop through active elements
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(iel,mdle,iflag,no,xsub,nrnodi,nodesi,nodesl,nvar,ndof, &
!$OMP            nr_elem_nodes,nface_orient,iv,ie,ifc,nod,nedge_orient, &
!$OMP            ntype,norder,zdofH,zdofE,zdofQ,tempH,tempE,tempV,ind)  &
!$OMP    SCHEDULE(DYNAMIC) REDUCTION(+:nr_up_elem)
      do iel=1,NRELES_SUBD
!
         mdle = ELEM_SUBD(iel)
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
            if (NODES(nod)%visit.eq.0 .and. is_Dirichlet(nod)) then
               nod_flg = .true.
               goto 100
            endif
         enddo
!
         nodesl(1:27) = nodes_elem(1:27,iel)
!
!     ...update the number of processed elements
         nr_up_elem = nr_up_elem+1
!
!-----------------------------------------------------------------------
!     ...Step 1: Update   V E R T   dof for Dirichlet nodes
         do iv=1,nvert(ntype)
            nod = abs(nodesl(iv))
            if (.not.associated(NODES(nod)%dof))       cycle
            if (.not.associated(NODES(nod)%dof%zdofH)) cycle
            if (NODES(nod)%visit.eq.1) cycle
            if (is_Dirichlet_attr(nod,CONTIN)) then
#if DEBUG_MODE
               if (iprint.eq.1) write(*,7010) mdle,iv,nod
          7010 format('update_Ddof: CALLING dhpvert FOR mdle,iv,nod = ',i8,i2,i8)
#endif
!
               nvar = size(NODES(nod)%dof%zdofH,1)
               ndof = size(NODES(nod)%dof%zdofH,2)
!
               call dhpvert(mdle,iflag,no,xsub(1:3,iv),NODES(nod)%case, &
                            NODES(nod)%bcond, tempH(1:nvar,1:ndof))
!
!$OMP CRITICAL
               if (NODES(nod)%visit.eq.0) then
                  NODES(nod)%dof%zdofH(:,:,N_COMS) = tempH(1:nvar,1:ndof)
                  NODES(nod)%visit = 1
               endif
!$OMP END CRITICAL
!
            endif
         enddo
!
!-----------------------------------------------------------------------
!
!     ...Step 2: Update   E D G E   dof for Dirichlet nodes
         call find_orient(mdle, nedge_orient,nface_orient)
         call find_order(mdle, norder)
!
!     ...compute solution dofs (need for H1 update)
         call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!
         do ie=1,nedge(ntype)
            ind = nvert(ntype)+ie
            nod = abs(nodesl(ind))
!
            if (.not.associated(NODES(nod)%dof)) cycle
            if (NODES(nod)%visit.ne.0)     cycle
            if (is_Dirichlet_attr(nod,CONTIN)) then
!           ...update H1 Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofH)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7020) mdle,ie,nod
             7020 format('update_Ddof: CALLING dhpedgeH FOR mdle,ie,nod = ',i8,i2,i8)
#endif
!
                  nvar = size(NODES(nod)%dof%zdofH,1)
                  ndof = size(NODES(nod)%dof%zdofH,2)
!
                  call dhpedgeH(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ie,    &
                                zdofH, tempH(1:nvar,1:ndof))
!
!$OMP CRITICAL
                  if (NODES(nod)%visit.eq.0) then
                     NODES(nod)%dof%zdofH(:,:,N_COMS) = tempH(1:nvar,1:ndof)
                     NODES(nod)%visit = 1
                  endif
!$OMP END CRITICAL
!
               else
                  NODES(nod)%visit=1
               endif
!
            endif
            if (is_Dirichlet_attr(nod,TANGEN)) then
!           ...update H(curl) Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofE)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7030) mdle,ie,nod
             7030 format('update_Ddof: CALLING dhpedgeE FOR mdle,ie,nod = ',i8,i2,i8)
#endif
!
                  nvar = size(NODES(nod)%dof%zdofE,1)
                  ndof = size(NODES(nod)%dof%zdofE,2)
!
                  call dhpedgeE(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ie,    &
                                tempE(1:nvar,1:ndof))
!
!$OMP CRITICAL
                  if (NODES(nod)%visit.lt.2) then
                     NODES(nod)%dof%zdofE(:,:,N_COMS) = tempE(1:nvar,1:ndof)
                     NODES(nod)%visit = 2
                  endif
!$OMP END CRITICAL
!
               else
                  NODES(nod)%visit=1
               endif
            endif
         enddo
!
!-----------------------------------------------------------------------
!  ...Step 3: Update   F A C E   dof for Dirichlet nodes
!
!     ...compute solution dofs (needed for H1 and Hcurl update)
         call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!
!     ...loop over faces
         do ifc=1,nface(ntype)
!        ...get local node number
            ind = nvert(ntype)+nedge(ntype)+ifc
!        ...get global node number
            nod = abs(nodesl(ind))
            if (.not.associated(NODES(nod)%dof)) cycle
            if (NODES(nod)%visit.ne.0)     cycle
            if (is_Dirichlet_attr(nod,CONTIN)) then
!           ...update H1 Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofH)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7040) mdle,ifc,nod
             7040 format('update_Ddof: CALLING dhpfaceH FOR mdle,ifc,nod = ',i8,i2,i8)
#endif
!
                  nvar = size(NODES(nod)%dof%zdofH,1)
                  ndof = size(NODES(nod)%dof%zdofH,2)
!
                  call dhpfaceH(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ifc,   &
                                zdofH, tempH(1:nvar,1:ndof))
!
!$OMP CRITICAL
                  if (NODES(nod)%visit.eq.0) then
                     NODES(nod)%dof%zdofH(:,:,N_COMS) = tempH(1:nvar,1:ndof)
                     NODES(nod)%visit = 1
                  endif
!$OMP END CRITICAL
!
               else
                  NODES(nod)%visit=1
               endif
            endif
            if (is_Dirichlet_attr(nod,TANGEN)) then
!           ...update H(curl) Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofE)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7050) mdle,ifc,nod
             7050 format('update_Ddof: CALLING dhpfaceE FOR mdle,ifc,nod = ',i8,i2,i8)
#endif
!
                  nvar = size(NODES(nod)%dof%zdofE,1)
                  ndof = size(NODES(nod)%dof%zdofE,2)
!
                  call dhpfaceE(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ifc,   &
                                zdofE, tempE(1:nvar,1:ndof))
!
!$OMP CRITICAL
                  if (NODES(nod)%visit.lt.2) then
                     NODES(nod)%dof%zdofE(:,:,N_COMS) = tempE(1:nvar,1:ndof)
                     NODES(nod)%visit = 2
                  endif
!$OMP END CRITICAL
!
               else
                  NODES(nod)%visit=1
               endif
            endif
            if (is_Dirichlet_attr(nod,NORMAL)) then
!           ...update H(div) Dirichlet dofs
               if (associated(NODES(nod)%dof%zdofV)) then
#if DEBUG_MODE
                  if (iprint.eq.1) write(*,7060) mdle,ifc,nod
             7060 format('update_Ddof: CALLING dhpfaceV FOR mdle,ifc,nod = ',i8,i2,i8)
#endif
!
                  nvar = size(NODES(nod)%dof%zdofV,1)
                  ndof = size(NODES(nod)%dof%zdofV,2)
!
                  call dhpfaceV(mdle,iflag,no,xsub,                     &
                                ntype,NODES(nod)%case,NODES(nod)%bcond, &
                                nedge_orient,nface_orient,norder,ifc,   &
                                tempV(1:nvar,1:ndof))
!
!$OMP CRITICAL
                  if (NODES(nod)%visit.lt.3) then
                     NODES(nod)%dof%zdofV(:,:,N_COMS) = tempV(1:nvar,1:ndof)
                     NODES(nod)%visit = 3
                  endif
!$OMP END CRITICAL
!
               else
                  NODES(nod)%visit=1
               endif
            endif
!     ...end of loop over faces
         enddo
!
         NODES(mdle)%visit=1
!
!  ...end of loop through elements
 100  continue
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
         write(*,*) 'update_Ddof: mesh is not distributed (or host mesh), but flag = ',nod_flg
         stop
      endif
      goto 70
   endif
!
   !write(*,4010) RANK,nod_flg
 4010 format('update_Ddof: [',I4,'], nod_flag=',L2)
!
!..skip creating node list if no DOFs are needed from other subdomains
   j_loc = 0
   if (.not. nod_flg) goto 60
!
!..fill local node list
   loc_max = 1000; allocate(nod_loc(loc_max))
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      ntype = NODES(mdle)%ntype
!
!  ...skip if the element has already been processed
      if (NODES(mdle)%visit.ge.1) cycle
!
!  ...check if irregular nodes have been updated; add to list if not
      nrnodi = nrnodes_irreg(iel)
      nodesi(1:nrnodi) = nodes_irreg(1:nrnodi,iel)
      do i=1,nrnodi
         nod = nodesi(i)
!
!     ...check if the node has been updated
         if ( NODES(nod)%visit.eq.0 .and. is_Dirichlet(nod) ) then
            NODES(nod)%visit = -1
            if(j_loc .ge. loc_max) then
               loc_max = loc_max*2
               allocate(nod_tmp(loc_max))
               nod_tmp(1:j_loc) = nod_loc(1:j_loc)
               call move_alloc(nod_tmp, nod_loc)
            endif
!
            j_loc = j_loc+1
            nod_loc(j_loc) = nod
!
         endif
      enddo
   enddo
!
!..consistency check
   if (j_loc .eq. 0) then
      write(*,*) 'update_Ddof: j_loc = 0, but flag = ',nod_flg
      stop
   endif
!
   !write(*,5010) RANK,j_loc
   5010 format('update_Ddof: [',I4,'], jloc=',I6)
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
      if (NODES(nod)%visit.ge.1) nod_rnk(i) = RANK
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
      tag = mod(nod_glb(i),200000)
      do while (i > k)
         rcv = rcv+1
         k = k+nod_cnt(rcv+1)
      enddo
      nod = nod_glb(i)
      call find_ndof(nod, ndofH,ndofE,ndofV,ndofQ)
      icase = NODES(Nod)%case
      nvarH = NREQNH(icase)*NRRHS
      nvarE = NREQNE(icase)*NRRHS
      nvarV = NREQNV(icase)*NRRHS
!  ...check whether supplier
      if (src .eq. RANK) then
         !write(*,6010) RANK,'SENDING  ',nod,S_Type(NODES(nod)%ntype)
         if (nvarH > 0 .and. ndofH > 0) then
            count = ndofH*nvarH; buf => NODES(nod)%dof%zdofH(:,:,N_COMS)
            call MPI_SEND(buf,count,MPI_VTYPE,rcv,tag,MPI_COMM_WORLD,ierr)
         endif
         if (nvarE > 0 .and. ndofE > 0) then
            count = ndofE*nvarE; buf => NODES(nod)%dof%zdofE(:,:,N_COMS)
            call MPI_SEND(buf,count,MPI_VTYPE,rcv,tag,MPI_COMM_WORLD,ierr)
         endif
         if (nvarV > 0 .and. ndofV > 0) then
            count = ndofV*nvarV; buf => NODES(nod)%dof%zdofV(:,:,N_COMS)
            call MPI_SEND(buf,count,MPI_VTYPE,rcv,tag,MPI_COMM_WORLD,ierr)
         endif
!  ...check whether receiver
      elseif (rcv .eq. RANK) then
         !write(*,6010) RANK,'RECEIVING',nod,S_Type(NODES(nod)%ntype)
         if (nvarH > 0 .and. ndofH > 0) then
            count = ndofH*nvarH; buf => NODES(nod)%dof%zdofH(:,:,N_COMS)
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         endif
         if (nvarE > 0 .and. ndofE > 0) then
            count = ndofE*nvarE; buf => NODES(nod)%dof%zdofE(:,:,N_COMS)
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         endif
         if (nvarV > 0 .and. ndofV > 0) then
            count = ndofV*nvarV; buf => NODES(nod)%dof%zdofV(:,:,N_COMS)
            call MPI_RECV(buf,count,MPI_VTYPE,src,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         endif
         NODES(nod)%visit = 3
      endif
 6010 format('update_Ddof: [',I4,'] ',A,' nod=',I6,', type=',A)
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
!..wait for non-blocking send/recv to finish (use MPI_ISEND/MPI_IRECV)
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
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,8010) end_time-start_time
 8010 format(' update_Ddof: ',f12.5,'  seconds',/)
   endif
!
end subroutine update_Ddof_omp
