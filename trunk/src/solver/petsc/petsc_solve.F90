!
#include "implicit_none.h"
#include <petsc/finclude/petscksp.h>
! -----------------------------------------------------------------------
!
!    routine name       - petsc_solve
!
! -----------------------------------------------------------------------
!
!    latest revision    - Feb 2020
!
!    purpose            - interface for distributed PETSc solvers
!                       - routine computes global stiffness matrix
!                         and global load vector and solves the global
!                         linear system with PETSc
!                       - the assembly is computed in parallel with one
!                         MPI process per subdomain, and OpenMP
!                         parallelization within each subdomain
!                       - this routine supports both computation with or
!                         without static condensation (uses module stc)
!
!    in                 - mtype:
!                                'P':  Positive definite (Herm./sym.)
!                                'H':  Hermitian (for complex)/
!                                      Symmetric (for real)
!                                      case for Lapack routines
!                                'G':  General case for Lapack routines
!
! -----------------------------------------------------------------------
subroutine petsc_solve(mtype)
!
   use data_structure3D, only: NRNODS, NRELES, NRELES_SUBD,       &
                               ELEM_ORDER, ELEM_SUBD,get_subd
   use assembly,         only: NR_RHS, MAXDOFM, MAXDOFS,          &
                               MAXbrickH, MAXmdlbH, NRHVAR,       &
                               MAXbrickE, MAXmdlbE, NREVAR,       &
                               MAXbrickV, MAXmdlbV, NRVVAR,       &
                               MAXbrickQ, NRQVAR,                 &
                               NEXTRACT, IDBC, ZDOFD, ZERO,       &
                               ALOC, BLOC, AAUX, ZAMOD, ZBMOD,    &
                               NR_PHYSA, MAXNODM
   use assembly_sc
   use control    ,  only: ISTC_FLAG
   use stc        , only:  HERM_STC,CLOC,                         &
                           stc_alloc,stc_dealloc,stc_get_nrdof
   use par_mesh   , only:  DISTRIBUTED,HOST_MESH
   use mpi_param  , only:  RANK,ROOT,NUM_PROCS
   use mpi_wrapper, only:  mpi_w_handle_err
   use MPI        , only:  MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,  &
                           MPI_INTEGER,MPI_INTEGER8,              &
                           MPI_REAL8,MPI_COMPLEX16,               &
                           MPI_COMM_WORLD,MPI_STATUS_SIZE
   use petscksp   , only:  VecCreateMPI,VecDuplicate,VecSet,      &
                           VecDestroy,MatCreate,MatDestroy,       &
                           VecSet,MatSetValues,MatSetSizes,       &
                           MatSetFromOptions,MatSetUp,            &
                           MatAssemblyBegin,MatAssemblyEnd,       &
                           VecAssemblyBegin,VecAssemblyEnd,       &
                           PETSC_VIEWER_STDOUT_WORLD,             &
                           PETSC_DECIDE,KSPSolve,KSPView,         &
                           MAT_FINAL_ASSEMBLY,ADD_VALUES,         &
                           INSERT_VALUES,PETSC_NULL_INTEGER,      &
                           VecScatterCreateToAll,VecScatterBegin, &
                           VecScatterEnd,VecScatterDestroy,       &
                           SCATTER_FORWARD,VecGetArrayReadF90,    &
                           VecRestoreArrayReadF90,tVec,           &
                           tVecScatter,KSPGetIterationNumber,     &
                           PetscPrintf,MatMPIAIJSetPreallocation, &
                           MAT_NEW_NONZERO_ALLOCATION_ERR,        &
                           PETSC_FALSE,VecGetOwnershipRange,      &
                           MatStashGetInfo,MatGetInfo,MAT_LOCAL,  &
                           MAT_INFO_SIZE,MAT_INFO_MALLOCS,        &
                           MAT_INFO_NZ_ALLOCATED,MAT_INFO_NZ_USED
   use petsc_w_ksp, only:  petsc_ksp_start,petsc_ksp_destroy,     &
                           petsc_ksp,petsc_A,petsc_rhs,petsc_sol
!
   implicit none
!
   character, intent(in) :: mtype
!
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!
!..number of local element dof for each physics variable
   integer, dimension(NR_PHYSA) :: nrdofi,nrdofb
! 
!..integer counters
   integer    :: nrdofm,nrdofc,nrnodm,nrdof,nrdof_mdl,ndof
   integer    :: iel,mdle,subd,idx,i,j,k,l,k1,k2,inod,nod,idec
!
!..MPI variables
   integer :: count,src,rcv,tag,ierr,nr_send,nr_recv
   integer :: reqs(2*NUM_PROCS),stats(MPI_STATUS_SIZE,2*NUM_PROCS)
!
!..PETSc variables
   PetscErrorCode petsc_ierr
   PetscScalar petsc_val ! can be real (8byte) or complex (16byte), depending on library linking
   PetscInt petsc_row, petsc_col, petsc_its
   PetscInt petsc_low, petsc_high, petsc_void
   PetscInt petsc_nstash, petsc_reallocs, petsc_nvoid
   PetscInt, allocatable :: petsc_dnz(:), petsc_onz(:)
   Vec petsc_sol_glb
   VecScatter petsc_ctx
   PetscScalar, pointer :: petsc_vec_ptr(:)
   real(8) :: petsc_info(MAT_INFO_SIZE)
   real(8) :: petsc_mallocs,petsc_nz_alloc,petsc_nz_used
   character*64 :: info_string
   character*8  :: fmt,val_string
!
!..Non-zero computation
   integer :: my_dnz, my_onz
   integer, allocatable :: NOD_INT(:,:)
   integer, allocatable :: NR_NOD_INT(:),NOD_DOF(:),NOD_LIST(:),NOD_VIS(:)
   integer :: NR_NOD_LIST
   integer, allocatable :: NOD_COMM(:,:)
   type ONZ_BUF
      integer, allocatable :: SEND_BUF(:)
      integer, allocatable :: RECV_BUF(:) 
   end type
   type(ONZ_BUF), allocatable :: ONZ(:)
   integer, allocatable :: TEMP_BUF(:)
!
!..dummy variables
   integer :: nvoid
   VTYPE   :: zvoid
! 
!..workspace for celem
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
!
!..64 bit non-zero entry counters
   integer(8) :: nnz,nnz_loc
   integer(8) :: elem_nnz_loc(NRELES_SUBD)
!
!..subdomain dof counters
   integer :: nrdof_subd(NUM_PROCS)
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time,time_stamp
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
   if ((.not. DISTRIBUTED) .or. HOST_MESH) then
      if (RANK .eq. ROOT) then
         write(*,*) 'petsc_solve: mesh is not distributed (or on host).'
         write(*,*) 'returning...'
      endif
      return
   endif
!
   select case(mtype)
      case('H')
         HERM_STC = .true.
      case default
         HERM_STC = .false.
   end select
!
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      write(*,1000)
1000  format(' petsc_solve: STARTED')
      write(*,*)
   endif
!
!..TODO multiple right-hand sides
   NR_RHS = 1
!
!..Initialize PETSc environment and allocate KSP data structure
   call petsc_ksp_start
!
! ----------------------------------------------------------------------
!  STEP 1 : 1ST LOOP THROUGH ELEMENTS, 1ST CALL TO CELEM TO GET INFO
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      if (RANK .eq. ROOT) write(*,1001)
 1001 format(' STEP 1 started : Get assembly info')
      start_time = MPI_Wtime()
   endif
!
!..allocate node ownership array
   allocate(NOD_OWN(NRNODS)); NOD_OWN = NUM_PROCS
!
!..compute local node ownership
!  (assumes node subdomains have previously been set)
!$OMP PARALLEL DO PRIVATE(nod,subd)
   do nod=1,NRNODS
      call get_subd(nod, subd)
      if (subd .eq. RANK) NOD_OWN(nod) = subd
   enddo
!$OMP END PARALLEL DO
!
!..compute global node ownership
   count = NRNODS
   call MPI_ALLREDUCE(MPI_IN_PLACE,NOD_OWN,count,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD, ierr)
   call mpi_w_handle_err(ierr,'MPI_ALLREDUCE')
!
   allocate(MAXDOFS(NR_PHYSA))
   MAXDOFS = 0; MAXDOFM = 0
!
   allocate(NFIRST_DOF(NRNODS)); NFIRST_DOF = -1
!
!..Step 1: determine the first dof offsets for active nodes
   nrdof = 0; nrdof_mdl = 0; idec = 1
!
!..non-zero counters for element offsets in distributed sparse stiffness matrix
!  using 64 bit integers nnz and nnz_loc
   nnz     = 0_8 ! global number of non-zeros in matrix (counts duplicate indices)
   nnz_loc = 0_8 ! local  number of non-zeros in matrix (counts duplicate indices)
   elem_nnz_loc(1:NRELES_SUBD) = 0_8 ! local element offsets for subdomain matrix
!
!..TODO: size NRNODS arrays are too expensive. We can use a NOD_VIS array to map to local NOD_LIST
   allocate(NOD_INT(250,NRNODS)); NOD_INT = 0
   allocate(NR_NOD_INT(NRNODS)); NR_NOD_INT = 0
   allocate(NOD_DOF(NRNODS)); NOD_DOF = 0
   allocate(NOD_LIST(NRNODS)); NOD_LIST = 0
   NR_NOD_LIST = 0
!..Every processor fills its row
!  each column specified the number of integer values that will be sent to that particular proc
   allocate(NOD_COMM(NUM_PROCS,NUM_PROCS)); NOD_COMM = 0
   allocate(ONZ(NUM_PROCS))
!
!..counter for temporary non-zero allocations that are assembled on another processor
   petsc_nstash = 0
!
!..compute offsets for owned nodes
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
!  ...get information from celem
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle,idec, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      else
         call celem(mdle,idec, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      endif
!
!  ...nrdofc = number of modified element dof after compression
!  ...k      = number of non-zero entries in element stiffness matrix
      k = nrdofc**2
!
!  ...subdomain counters for OpenMP
      elem_nnz_loc(iel) = nnz_loc
      nnz_loc = nnz_loc + int8(k)
!
!  ...for PETSc distributed sparse matrix assembly, we must compute the number of non-zero
!     entries per row (dof);
!
!  ...update the maximum number of local dof
      do i=1,NR_PHYSA
         MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
      enddo
!  ...update the maximum number of modified element dof in the expanded mode
      MAXDOFM = max0(MAXDOFM,nrdofm)
!
!  ...create list of node interaction
      do i=1,nrnodm
         nod = nodm(i)
         if (NOD_DOF(nod) .eq. 0) then
            NOD_DOF(nod) = ndofmH(i) + ndofmE(i) + ndofmV(i) + ndofmQ(i)
            if (NOD_DOF(nod) .eq. 0) cycle ! no interaction
            NR_NOD_LIST = NR_NOD_LIST + 1
            NOD_LIST(NR_NOD_LIST) = nod
         endif
         do j=1,nrnodm
            k = nodm(j)
            if (ndofmH(j) + ndofmE(j) + ndofmV(j) + ndofmQ(j) .eq. 0) cycle ! no interaction
            do l = 1,NR_NOD_INT(nod)
               if (NOD_INT(l,nod) .eq. k) then
                  k = 0
                  exit
               endif
            enddo
            if (k > 0) then
!           ...skip if node interaction with itself if owned by another subdomain
               if ((NOD_OWN(nod) .ne. RANK) .and. (nod .eq. k)) cycle
               NR_NOD_INT(nod) = NR_NOD_INT(nod) + 1
               if (NR_NOD_INT(nod) > 250) then
                  write(*,*) 'NR_NOD_INT(nod) > 250. stop.'; stop
               endif
               NOD_INT(NR_NOD_INT(nod),nod) = k
            endif
         enddo
      enddo
!
!  ...compute offsets for nodal dof
      do i=1,nrnodm
         nod = nodm(i)
!     ...avoid repetition within overlaps with other subdomains
         if (NOD_OWN(nod) .ne. RANK) cycle
!     ...avoid repetition within my subdomain
         if (NFIRST_DOF(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRST_DOF(nod) = nrdof
!     ...update the dof counter
         nrdof = nrdof + ndofmH(i) + ndofmE(i) + ndofmV(i)
      enddo
      if (.not. ISTC_FLAG) then
         nrdof = nrdof + ndofmQ(nrnodm)
      endif
!
!  ...compute number of bubble dof (nrdof_mdl)
      if (ISTC_FLAG) then
         call stc_get_nrdof(mdle, nrdofi,nrdofb)
         nrdof_mdl = nrdof_mdl + sum(nrdofb)
      endif
   enddo
!
!..compute subdomain offset
   nrdof_subd(1:NUM_PROCS) = 0
   nrdof_subd(RANK+1) = nrdof
   count = NUM_PROCS
   call MPI_ALLREDUCE(MPI_IN_PLACE,nrdof_subd,count,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD, ierr)
   call mpi_w_handle_err(ierr,'MPI_ALLREDUCE')
!
!..compute non-zero entries per row, accounting for own subdomain interaction
   allocate(petsc_dnz(nrdof_subd(RANK+1))); petsc_dnz = 0
   allocate(petsc_onz(nrdof_subd(RANK+1))); petsc_onz = 0
   do i=1,NR_NOD_LIST
      nod = NOD_LIST(i)
      if (NOD_DOF(nod) .eq. 0) then ! no interaction with any node
         write(*,*) 'NOD_DOF(nod) .eq. 0 !!!'; stop ! CHECK (should not happen)
      endif
!  ...if node is not owned, the node interactions need to be communicated to the owner
      if (NOD_OWN(nod) .ne. RANK) then
         j = NOD_OWN(nod) ! rank of the node owner
!     ...check current buffer size
         if (.not. allocated(ONZ(j+1)%SEND_BUF) ) allocate(ONZ(j+1)%SEND_BUF(10000))
         k = size(ONZ(j+1)%SEND_BUF)
         if (NOD_COMM(RANK+1,j+1)+(2+2*NR_NOD_INT(nod)) > k) then
            write(*,*) 'Reallocating SEND_BUF...'
            allocate(TEMP_BUF(2*k))
            TEMP_BUF(1:k) = ONZ(j+1)%SEND_BUF(1:k)
            call move_alloc(TEMP_BUF, ONZ(j+1)%SEND_BUF) 
         endif
!     ...fill buffer with this node's interaction with my own subdomain
!        data = (nod, NR_NOD_INT(nod), (nod_int,dof), (nod_int,dof), .. )
         ONZ(j+1)%SEND_BUF(NOD_COMM(RANK+1,j+1)+1) = nod
         !write(*,*) 'Adding node to send buffer: nod = ', nod
         k1 = 0
         do l = 1,NR_NOD_INT(nod)
            k = NOD_INT(l,nod)
            if (k .eq. j) then ! owner already knows about this node interaction
               write(*,*) 'k .eq. j !!!'; stop ! CHECK (should not happen)
            endif
            if (NOD_DOF(k) .eq. 0) then ! no interaction with this node
               write(*,*) 'NOD_DOF(k) .eq. 0 ...'; stop ! CHECK (should not happen)
            endif
            ONZ(j+1)%SEND_BUF(NOD_COMM(RANK+1,j+1)+2+k1+1) = k
            ONZ(j+1)%SEND_BUF(NOD_COMM(RANK+1,j+1)+2+k1+2) = NOD_DOF(k)
            petsc_nstash = petsc_nstash + NOD_DOF(k)*NOD_DOF(nod)
            k1 = k1+2
         enddo
         ONZ(j+1)%SEND_BUF(NOD_COMM(RANK+1,j+1)+2) = k1
!     ...increment counter in NOD_COMM list (by the amount of data to be sent)
         NOD_COMM(RANK+1,j+1) = NOD_COMM(RANK+1,j+1) + (2+k1)
         cycle
      endif
      my_dnz = 0; my_onz = 0
      do j = 1,NR_NOD_INT(nod)
         k = NOD_INT(j,nod)
         if (NOD_DOF(k) .eq. 0) then
            write(*,*) 'NOD_DOF(k) .eq. 0 !!!'; stop ! CHECK (should not happen)
         endif
         if (NOD_OWN(k) .eq. RANK) then
            my_dnz = my_dnz + NOD_DOF(k)
         else
            my_onz = my_onz + NOD_DOF(k)
         endif
      enddo
      if (NFIRST_DOF(nod) < 0) then
         write(*,*) 'NFIRST_DOF(nod) < 0. stop.'; stop ! CHECK (should not happen)
      endif
      if (NFIRST_DOF(nod)+NOD_DOF(nod) > nrdof_subd(RANK+1)) then
         write(*,*) 'NFIRST_DOF(nod)+NOD_DOF(nod) > nrdof_subd(RANK+1). stop.'; stop ! CHECK (should not happen)
      endif
      do j=1,NOD_DOF(nod)
         if (petsc_dnz(NFIRST_DOF(nod)+j) .ne. 0) then
            write(*,*) 'petsc_dnz(NFIRST_DOF(nod)+j) .ne. 0'; stop ! CHECK (should not happen)
         endif
      enddo
      petsc_dnz(NFIRST_DOF(nod)+1 : NFIRST_DOF(nod)+NOD_DOF(nod)) = my_dnz
      petsc_onz(NFIRST_DOF(nod)+1 : NFIRST_DOF(nod)+NOD_DOF(nod)) = my_onz
   enddo
!
!..communicate which off-diagonal contributions need to be exchanged
!  in order to account for non-zero entries from interaction between subdomains
   count = NUM_PROCS * NUM_PROCS
   call MPI_ALLREDUCE(MPI_IN_PLACE,NOD_COMM(1:NUM_PROCS,1:NUM_PROCS),count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
   call mpi_w_handle_err(ierr,'MPI_ALLREDUCE')
!
!..exchange off-diagonal contributions
   nr_send = 0; nr_recv = 0
   do i=1,NUM_PROCS
!  ...perform all send operations
      if (NOD_COMM(RANK+1,i) > 0) then
         count = NOD_COMM(RANK+1,i)
         rcv   = i-1 ! RANK of receiving proc
         tag   = count
         nr_send  = nr_send + 1
         call MPI_ISEND(ONZ(i)%SEND_BUF(1:count),count,MPI_INTEGER,rcv,tag,MPI_COMM_WORLD, reqs(nr_send),ierr)
         call mpi_w_handle_err(ierr,'MPI_ISEND')
      endif
   enddo
   do i=1,NUM_PROCS
!  ...perform all receive operations
      if (NOD_COMM(i,RANK+1) > 0) then
         count = NOD_COMM(i,RANK+1)
         src   = i-1 ! RANK of sending proc
         tag   = count
         nr_recv = nr_recv + 1
         allocate(ONZ(i)%RECV_BUF(count))
         call MPI_IRECV(ONZ(i)%RECV_BUF(1:count),count,MPI_INTEGER,src,tag,MPI_COMM_WORLD, reqs(nr_send+nr_recv),ierr)
         call mpi_w_handle_err(ierr,'MPI_IRECV')
      endif
   enddo
   !write(*,'(A,I2,A,I2)') '[RANK], nr_send = [',RANK,'],',nr_send
   !write(*,'(A,I2,A,I2)') '[RANK], nr_recv = [',RANK,'],',nr_recv
!..wait until all send and receive operations are completed
   call MPI_Waitall(nr_send+nr_recv,reqs(1:nr_send+nr_recv), stats(:,1:nr_send+nr_recv),ierr)
   call mpi_w_handle_err(ierr,'MPI_Waitall')
!
!..process the received node interactions
   do i=1,NUM_PROCS
      if (NOD_COMM(i,RANK+1) > 0) then
         count = NOD_COMM(i,RANK+1)
         l = 0
         do while (l < count)
            nod = ONZ(i)%RECV_BUF(l+1) ! node number
            !write(*,*) 'Processing node to from receive buffer: nod = ', nod
            k1  = ONZ(i)%RECV_BUF(l+2) ! number of received node interactions: (nod,ndof) pairs
            !write(*,*) 'Number of interactions: k1 = ',k1
            my_dnz = 0; my_onz = 0
            do j = 1,k1,2
               k    = ONZ(i)%RECV_BUF(l+2+j  ) ! fetch node interaction
               ndof = ONZ(i)%RECV_BUF(l+2+j+1) ! add number of dof of this interaction
               if (NOD_DOF(k) .eq. 0) NOD_DOF(k) = ndof
!           ...add node interaction if it has not yet been accounted for
               do k2 = 1,NR_NOD_INT(nod)
                  inod = NOD_INT(k2,nod)
                  if (k .eq. inod) then
                     k = 0; exit
                  endif 
               enddo
               if (k > 0) then
                  if (NOD_OWN(k) .eq. RANK) then
                     my_dnz = my_dnz + ndof
                  else
                     my_onz = my_onz + ndof
                  endif
                  NR_NOD_INT(nod) = NR_NOD_INT(nod) + 1
                  NOD_INT(NR_NOD_INT(nod),nod) = k
               endif
            enddo
            if (NFIRST_DOF(nod) < 0) then
               write(*,*) 'NFIRST_DOF(nod) < 0'; stop ! CHECK
            endif
            if (NOD_DOF(nod) .eq. 0) then
               write(*,*) 'NOD_DOF(nod) = 0'; stop ! CHECK
            endif
            petsc_dnz(NFIRST_DOF(nod)+1 : NFIRST_DOF(nod)+NOD_DOF(nod)) = &
            petsc_dnz(NFIRST_DOF(nod)+1 : NFIRST_DOF(nod)+NOD_DOF(nod)) + my_dnz
            petsc_onz(NFIRST_DOF(nod)+1 : NFIRST_DOF(nod)+NOD_DOF(nod)) = &
            petsc_onz(NFIRST_DOF(nod)+1 : NFIRST_DOF(nod)+NOD_DOF(nod)) + my_onz
            l = l + (2+k1)
         enddo
      endif
   enddo
!
!..deallocate auxiliary arrays
   do i=1,NUM_PROCS
      if (allocated(ONZ(i)%SEND_BUF)) deallocate(ONZ(i)%SEND_BUF)
      if (allocated(ONZ(i)%RECV_BUF)) deallocate(ONZ(i)%RECV_BUF)
   enddo
   deallocate(NOD_LIST,NOD_DOF,NOD_INT,NR_NOD_INT,NOD_COMM,ONZ)
   write(*,'(A,I2,A,I8)') '[RANK], sum(petsc_dnz) = [',RANK,'],',sum(petsc_dnz)
   write(*,'(A,I2,A,I8)') '[RANK], sum(petsc_onz) = [',RANK,'],',sum(petsc_onz)
   !write(*,'(A,I2,A,I8)') '[RANK], sum(dnz + onz) = [',RANK,'],',sum(petsc_dnz) + sum(petsc_onz)
!
!..calculate prefix sum for global offsets
   nrdof = 0
   do i = 1,RANK
      nrdof = nrdof + nrdof_subd(i)
   enddo
!$OMP PARALLEL DO
   do i = 1,NRNODS
      if (NOD_OWN(i) .eq. RANK) NFIRST_DOF(i) = NFIRST_DOF(i) + nrdof
   enddo
!$OMP END PARALLEL DO
!
!..communicate offsets (to receive offsets for non-owned nodes within subdomain)
   count = NRNODS
   call MPI_ALLREDUCE(MPI_IN_PLACE,NFIRST_DOF,count,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD, ierr)
   call mpi_w_handle_err(ierr,'MPI_ALLREDUCE')
!
!..calculate total number of (interface) dofs
   nrdof = 0
   do i = 1,NUM_PROCS
      nrdof = nrdof + nrdof_subd(i)
   enddo
!
!..compute total number of condensed bubble dofs
   count = 1
   call MPI_ALLREDUCE(MPI_IN_PLACE,nrdof_mdl,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
   call mpi_w_handle_err(ierr,'MPI_ALLREDUCE')
!
!..compute total number of non-zeros in global matrix
   count = 1
   call MPI_ALLREDUCE(nnz_loc,nnz,count,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD, ierr)
   call mpi_w_handle_err(ierr,'MPI_ALLREDUCE')
!
!..total number of (interface) dof is nrdof
   NRDOF_CON = nrdof
   NRDOF_TOT = nrdof + nrdof_mdl
!
   deallocate(NOD_OWN)
!
   if (nrdof .eq. 0) then
      deallocate(MAXDOFS,NFIRST_DOF)
      if (RANK .eq. ROOT) write(*,*) 'petsc_solve: nrdof = 0. returning.'
      return
   endif
!
! ----------------------------------------------------------------------
!  END OF STEP 1
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      end_time = MPI_Wtime()
      Mtime(1) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1002) Mtime(1)
 1002 format(' STEP 1 finished: ',f12.5,'  seconds',/)
   endif 
!
! ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM 
! ----------------------------------------------------------------------
!
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call mpi_w_handle_err(ierr,'MPI_BARRIER')
   if (RANK .eq. ROOT) then
      write(*,2010) '[', RANK, '] Number of dof  : nrdof_con = ', NRDOF_CON
      write(*,2010) '[', RANK, ']                  nrdof_tot = ', NRDOF_TOT
      write(*,2010) '[', RANK, '] Total non-zeros: nnz       = ', nnz
   endif
   write(*,2010) '[', RANK, '] Local non-zeros: nnz_loc   = ', nnz_loc
2010 format(A,I4,A,I12)
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      if (RANK .eq. ROOT) write(*,1003)
 1003 format(/,' STEP 2 started : Global Assembly')
      start_time = MPI_Wtime()
   endif
!
!..memory allocation for PETSc solver
!..create distributed load vector (right-hand side)
   call VecCreateMPI(MPI_COMM_WORLD,nrdof_subd(RANK+1),nrdof, petsc_rhs,petsc_ierr); CHKERRQ(petsc_ierr)
!..create distributed solution vector
   call VecDuplicate(petsc_rhs, petsc_sol,petsc_ierr); CHKERRQ(petsc_ierr)
!
   if (RANK .eq. ROOT) then
      do i=1,NUM_PROCS
         write(*,5678) i-1,nrdof_subd(i)
         5678 format('[',I2,']: nrdof_subd = ',I7)
      enddo
   endif
!
   call VecGetOwnershipRange(petsc_rhs, petsc_low,petsc_high,petsc_ierr); CHKERRQ(petsc_ierr)
   write(*,5679) RANK,petsc_low,petsc_high
 5679 format('[',I2,']: ',I8,' to ',I8)
!
!..create distributed sparse stiffness matrix
   call MatCreate(MPI_COMM_WORLD, petsc_A,petsc_ierr); CHKERRQ(petsc_ierr)
   call MatSetSizes(petsc_A,nrdof_subd(RANK+1),nrdof_subd(RANK+1),nrdof,nrdof, petsc_ierr); CHKERRQ(petsc_ierr)
   call MatSetFromOptions(petsc_A, petsc_ierr); CHKERRQ(petsc_ierr)
   !call MatSetOption(petsc_A,MAT_SYMMETRIC,PETSC_TRUE, petsc_ierr); CHKERRQ(petsc_ierr);
   !call MatSetOption(petsc_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE, petsc_ierr); CHKERRQ(petsc_ierr)
!..preallocate sparse matrix
   !petsc_dnz(1:nrdof_subd(RANK+1)) = nrdof_subd(RANK+1) ! preallocate with upper bound
   !petsc_onz(1:nrdof_subd(RANK+1)) = nrdof-nrdof_subd(RANK+1) ! preallocate with upper bound
   petsc_void = 1
   call MatMPIAIJSetPreallocation(petsc_A,petsc_void,petsc_dnz,petsc_void,petsc_onz, petsc_ierr); CHKERRQ(petsc_ierr)
   if (petsc_nstash > 0) petsc_nstash = MAX(2*petsc_nstash, 10000)
   call MatStashSetInitialSize(petsc_A,petsc_nstash,petsc_void, petsc_ierr); CHKERRQ(petsc_ierr)
   deallocate(petsc_dnz,petsc_onz)
!
!..memory allocation for static condensation
   call stc_alloc
!
!..assemble global stiffness matrix
!..loop through elements
   idec = 2
!
!$OMP PARALLEL                                  &
!$OMP PRIVATE(nrdofs,nrdofm,nrdofc,nodm,nrnodm, &
!$OMP         ndofmH,ndofmE,ndofmV,ndofmQ,      &
!$OMP         i,j,k,k1,k2,l,nod,ndof,mdle,subd, &
!$OMP         petsc_row,petsc_col,petsc_val)
   allocate(NEXTRACT(MAXDOFM))
   allocate(IDBC(MAXDOFM))
   allocate(ZDOFD(MAXDOFM,NR_RHS))
   allocate(BLOC(NR_PHYSA))
   allocate(AAUX(NR_PHYSA))
   allocate(ALOC(NR_PHYSA,NR_PHYSA))
   do i=1,NR_PHYSA
      BLOC(i)%nrow = MAXDOFS(i)
      BLOC(i)%ncol = NR_RHS
      allocate(BLOC(i)%array(MAXDOFS(i),NR_RHS))
      do j=1,NR_PHYSA
         ALOC(i,j)%nrow = MAXDOFS(i)
         ALOC(i,j)%ncol = MAXDOFS(j)
         allocate(ALOC(i,j)%array(MAXDOFS(i),MAXDOFS(j)))
      enddo
      AAUX(i)%nrow = MAXDOFM
      AAUX(i)%ncol = MAXDOFS(i)
      allocate(AAUX(i)%array(MAXDOFM,MAXDOFS(i)))
   enddo
   allocate(ZBMOD(MAXDOFM,NR_RHS))
   allocate(ZAMOD(MAXDOFM,MAXDOFM))
   allocate(LCON(MAXDOFM))
   allocate(ZLOAD(MAXDOFM))
   allocate(ZTEMP(MAXDOFM**2))
!
!$OMP DO                 &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle,idec, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      else
         call celem(mdle,idec, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      endif
!
!  ...determine local to global dof connectivities
      l=0 ! element dof counter
!
!  ...H1 dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmH(i)
            l=l+1
            LCON(l) = NFIRST_DOF(nod)+j
         enddo
      enddo
!  ...H(curl) dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmE(i)
            l=l+1
            LCON(l) = NFIRST_DOF(nod)+ndofmH(i)+j
         enddo
      enddo
!  ...H(div) dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmV(i)
            l=l+1
            LCON(l) = NFIRST_DOF(nod)+ndofmH(i)+ndofmE(i)+j
         enddo
      enddo
!  ...L2 dof
      if (.not. ISTC_FLAG) then
         nod = nodm(nrnodm)
         do j=1,ndofmQ(nrnodm)
            l=l+1
            LCON(l) = NFIRST_DOF(nod)+ndofmH(nrnodm)+ndofmE(nrnodm)+ndofmV(nrnodm)+j
         enddo
      endif
!
!  ...number of element (interface) dof
      ndof = l
!
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
!  ...omp critical is needed for vecsetvalue,matsetvalue
!  ...TODO: assemble one block with OpenMP, then call setvalues in critical section
!$OMP CRITICAL
      do k1=1,ndof
!     ...global dof is:
         petsc_row = LCON(k1)-1 ! petsc indices start with 0
!     ...Assemble global load vector
         petsc_val = ZLOAD(k1)
         call VecSetValue(petsc_rhs,petsc_row,petsc_val,ADD_VALUES, petsc_ierr)!; CHKERRQ(petsc_ierr)
      enddo
      do k1=1,ndof
!     ...global dof is:
         petsc_row = LCON(k1)-1 ! petsc indices start with 0
!     ...loop through dof `to the right'
         do k2=1,ndof
!        ...global dof is:
            petsc_col = LCON(k2)-1 ! petsc indices start with 0
!        ...assemble
!        ...note: repeated indices are summed automatically by PETSc
            !elem_nnz_loc(iel) = elem_nnz_loc(iel) + 1_8
            k = (k1-1)*ndof + k2
            petsc_val = ZTEMP(k)
            call MatSetValue(petsc_A,petsc_row,petsc_col,petsc_val,ADD_VALUES, petsc_ierr)!; CHKERRQ(petsc_ierr)
         enddo
      enddo
!$OMP END CRITICAL
!
      CLOC(iel)%ni = ndof
      allocate(CLOC(iel)%con(ndof))
      CLOC(iel)%con = LCON(1:ndof)
!..end of loop through elements
   enddo
!$OMP END DO
!
   do i=1,NR_PHYSA
      deallocate(BLOC(i)%array)
      do j=1,NR_PHYSA
         deallocate(ALOC(i,j)%array)
      enddo
         deallocate(AAUX(i)%array)
   enddo
!
   deallocate(NEXTRACT,IDBC,ZDOFD,BLOC,AAUX,ALOC)
   deallocate(ZBMOD,ZAMOD,LCON,ZLOAD,ZTEMP)
!$OMP END PARALLEL
!
   deallocate(MAXDOFS,NFIRST_DOF)
!
!..PETSc Load VectorAssembly
   if (RANK .eq. ROOT) write(*,*) ' PETSc Assembly...'
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      time_stamp = MPI_Wtime()
   endif
   call VecAssemblyBegin(petsc_rhs, petsc_ierr); CHKERRQ(petsc_ierr)
   call VecAssemblyEnd  (petsc_rhs, petsc_ierr); CHKERRQ(petsc_ierr)
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      time_stamp = MPI_Wtime()-time_stamp
      if (RANK .eq. ROOT) write(*,3001) time_stamp
 3001 format(' - VecAssembly : ',f12.5,'  seconds')
   endif
!..PETSc Stiffness Matrix Assembly
   write(*,'(A,I2,A,I8)') '[RANK], petsc_nstash   = [',RANK,'],',petsc_nstash
   call MatStashGetInfo(petsc_A, petsc_nstash,petsc_reallocs,petsc_void,petsc_nvoid,petsc_ierr); CHKERRQ(petsc_ierr)
   write(*,2348) RANK,petsc_nstash,petsc_reallocs
   2348 format('[',I2,']',': nstash = ',I6,', reallocs = ',I6)
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      time_stamp = MPI_Wtime()
   endif
   call MatAssemblyBegin(petsc_A,MAT_FINAL_ASSEMBLY, petsc_ierr); CHKERRQ(petsc_ierr)
   call MatAssemblyEnd  (petsc_A,MAT_FINAL_ASSEMBLY, petsc_ierr); CHKERRQ(petsc_ierr)
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      time_stamp = MPI_Wtime()-time_stamp
      if (RANK .eq. ROOT) write(*,3002) time_stamp
 3002 format(' - MatAssembly : ',f12.5,'  seconds')
   endif
   call MatGetInfo(petsc_A,MAT_LOCAL, petsc_info,petsc_ierr); CHKERRQ(petsc_ierr)
   petsc_mallocs  = petsc_info(MAT_INFO_MALLOCS)
   petsc_nz_alloc = petsc_info(MAT_INFO_NZ_ALLOCATED)
   petsc_nz_used  = petsc_info(MAT_INFO_NZ_USED)
   write(*,2349) RANK,'petsc_mallocs  = ',INT(petsc_mallocs)
   write(*,2349) RANK,'petsc_nz_alloc = ',INT(petsc_nz_alloc)
   write(*,2349) RANK,'petsc_nz_used  = ',INT(petsc_nz_used)
   2349 format('[',I2,']',': ',A,I8)
   if (INT(petsc_nz_used) .ne. INT(petsc_nz_alloc)) then
      write(*,2349) RANK,'Preallocation was not accurate: diff = ', &
                    (INT(petsc_nz_alloc)-INT(petsc_nz_used))
   endif
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      end_time = MPI_Wtime()
      Mtime(2) =  end_time-start_time
      if (RANK .eq. ROOT) write(*,1004) Mtime(2)
 1004 format(' STEP 2 finished: ',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  STEP 3: call mumps to solve the linear system
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      if (RANK .eq. ROOT) write(*,1009)
 1009 format(' STEP 3 started : Solve')
      start_time = MPI_Wtime()
   endif
!
!..
   call KSPSetOperators(petsc_ksp,petsc_A,petsc_A, petsc_ierr); CHKERRQ(petsc_ierr)
!
!..
   !call KSPSetFromOptions(...)
   !call KSPSetType(petsc_ksp, KSPType petsc_method)
   ! e.g., KSPCG, KSPGMRES, ...
   ! default for CG is Hermitian positive definite
   !call KSPCGSetType(petsc_ksp, KSPCGType KSP_CG_SYMMETRIC)
!
   select case(mtype)
      case('P')
         ! CG
      case default
         ! GMRES
   end select
!
!..default solver is restarted GMRES
   call KSPSolve(petsc_ksp,petsc_rhs, petsc_sol,petsc_ierr); CHKERRQ(petsc_ierr)
!
!..print number of iterations
   call KSPGetIterationNumber(petsc_ksp, petsc_its,petsc_ierr); CHKERRQ(petsc_ierr);
   fmt = '(I5)'
   write (val_string,fmt) petsc_its
   info_string='KSPSolve: number of iterations = '//trim(val_string)
   call PetscPrintf(MPI_COMM_WORLD,info_string, petsc_ierr); CHKERRQ(petsc_ierr);
   if (RANK .eq. ROOT) write(*,*)
   if (RANK .eq. ROOT) write(*,*)
!
!..print solver info
   !call KSPView(petsc_ksp,PETSC_VIEWER_STDOUT_WORLD, petsc_ierr); CHKERRQ(petsc_ierr)
!
! ----------------------------------------------------------------------
!  END OF STEP 3
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      end_time = MPI_Wtime()
      Mtime(3) =  end_time-start_time
      if (RANK .eq. ROOT) write(*,1010) Mtime(3)
 1010 format(' STEP 3 finished: ',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 4 : STORE SOLUTION IN THE DATASTRACTURE
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      if (RANK .eq. ROOT) write(*,1011)
 1011 format(' STEP 4 started : Store the solution')
      start_time = MPI_Wtime()
   endif
!
!..broadcast global solution from host to other processes
   call VecScatterCreateToAll(petsc_sol, petsc_ctx,petsc_sol_glb,petsc_ierr); CHKERRQ(petsc_ierr)
   call VecScatterBegin(petsc_ctx,petsc_sol,petsc_sol_glb,INSERT_VALUES,SCATTER_FORWARD, petsc_ierr); CHKERRQ(petsc_ierr)
   call VecScatterEnd  (petsc_ctx,petsc_sol,petsc_sol_glb,INSERT_VALUES,SCATTER_FORWARD, petsc_ierr); CHKERRQ(petsc_ierr)
   call VecScatterDestroy(petsc_ctx, petsc_ierr); CHKERRQ(petsc_ierr)
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      time_stamp = MPI_Wtime()-start_time
      if (RANK .eq. ROOT) write(*,3004) time_stamp
 3004 format(' - Broadcast: ',f12.5,'  seconds')
   endif
!
   call VecGetArrayReadF90(petsc_sol_glb, petsc_vec_ptr,petsc_ierr); CHKERRQ(petsc_ierr)
!
   ndof = 0
!$OMP PARALLEL
!$OMP DO REDUCTION(MAX:ndof)
   do iel=1,NRELES_SUBD
      if (CLOC(iel)%ni > ndof) ndof = CLOC(iel)%ni
   enddo
!$OMP END DO
   allocate(ZSOL_LOC(ndof))
!$OMP DO PRIVATE(i,k1) SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      ZSOL_LOC=ZERO
      do k1=1,CLOC(iel)%ni
         i = CLOC(iel)%con(k1)
         ZSOL_LOC(k1) = petsc_vec_ptr(i)
      enddo
      deallocate(CLOC(iel)%con)
      call solout(iel,ndof,NR_RHS,ZERO,ZSOL_LOC)
   enddo
!$OMP END DO
   deallocate(ZSOL_LOC)
!$OMP END PARALLEL
!
   call VecRestoreArrayReadF90(petsc_sol_glb, petsc_vec_ptr,petsc_ierr); CHKERRQ(petsc_ierr)
   call VecDestroy(petsc_sol_glb,petsc_ierr); CHKERRQ(petsc_ierr)
!
! ----------------------------------------------------------------------
!  END OF STEP 4
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call mpi_w_handle_err(ierr,'MPI_BARRIER')
      end_time = MPI_Wtime()
      Mtime(4) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1012) Mtime(4)
 1012 format(' STEP 4 finished: ',f12.5,'  seconds',/)
   endif
!
   call stc_dealloc
!
   call VecDestroy(petsc_rhs,petsc_ierr); CHKERRQ(petsc_ierr)
   call VecDestroy(petsc_sol,petsc_ierr); CHKERRQ(petsc_ierr)
   call MatDestroy(petsc_A  ,petsc_ierr); CHKERRQ(petsc_ierr)
!
!..Finalize PETSc environment and destroy the KSP instance
   call petsc_ksp_destroy
!
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call mpi_w_handle_err(ierr,'MPI_BARRIER')
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .ge. 1)) then
      write(*,1013) sum(Mtime(1:4))
 1013 format(' petsc_solve FINISHED: ',f12.5,'  seconds',/)
   endif
!
!
end subroutine petsc_solve
