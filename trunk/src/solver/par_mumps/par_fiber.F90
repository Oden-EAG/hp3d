!
#include "typedefs.h"
!
#if HP3D_USE_INTEL_MKL
!
! -----------------------------------------------------------------------
!
!    routine name       - par_fiber
!
! -----------------------------------------------------------------------
!
!    latest revision    - Mar 2020
!
!    purpose            - interface for distributed MUMPS solver
!                       - routine solves distributed coupled interface
!                         problem specified in the input MUMPS instance
!                       - this solver splits the problem into several
!                         subproblems of specified size and solves each
!                         one separately (essentially performing static
!                         condensation onto the separator dofs)
!                       - the routine can be applied recursively in a
!                         divide-and-conquer fashion leading to a
!                         recursive nested dissection solve
!
!    in                 - mumps: MUMPS instance (see par_mumps.F90)
!                       - nrdof: interface dofs owned by each processor
!                                for the 'mumps' input instance
!                                (array of size 'nproc')
!                       - nproc: number of processors participating in
!                                the distributed solve of 'mumps' input
!                       - level: recursion counter (nested level)
!
!               note: this routine makes several assumptions about the
!                     structure (partitioning) of the input problem;
!                     it is only meant to perform recursive nested
!                     dissection on the waveguide problems that have
!                     been partitioned with the custom fiber partitioner
!                     and where subdomain interior dofs have been
!                     eliminated previously (via par_nested subroutine)
!
! -----------------------------------------------------------------------
recursive subroutine par_fiber(mumps,nrdof,nproc,level)
!
   use assembly  , only: NR_RHS
   use MPI       , only: MPI_UNDEFINED,MPI_STATUS_IGNORE,      &
                         MPI_REAL8,MPI_COMPLEX16,MPI_INTEGER,  &
                         MPI_IN_PLACE,MPI_SUM,MPI_Wtime
   use mpi_param , only: RANK,ROOT
   use parameters, only: ZERO,ZONE
   use stc       , only: stc_bwd
   use par_mumps
   use mkl_spblas
!
   implicit none
!
#if C_MODE
   type (ZMUMPS_STRUC), intent(inout) :: mumps
   type (ZMUMPS_STRUC) :: mumps_sub, mumps_int
#else
   type (DMUMPS_STRUC), intent(inout) :: mumps
   type (DMUMPS_STRUC) :: mumps_sub, mumps_int
#endif
!
!..number of owned dofs per MPI proc
!  (owned from the viewpoint of the nested dissection solver)
   integer, intent(in) :: nproc
   integer, intent(in) :: nrdof(nproc)
   integer, intent(in) :: level
!
!..set parameter to indicate the size of each subproblem
   integer, parameter :: mSUB_PROCS = 4
!
!..set parameter to indicate whether more than one level of
!  nested dissection should be performed
   logical, parameter :: RECURSIVE_NESTED = .true.
!
!..aux variables for subcommunicators
   integer :: mRANK,mPROCS,mSUB_RANK,mINT_PROCS,mINT_RANK
   integer :: mpi_comm_sub,mpi_comm_int
   integer :: group,key,color,count,src,rcv,tag,ierr
   integer :: nrdof_int(nproc/mSUB_PROCS)
!
!..aux variables for subproblem assembly
   integer :: dof_off,dof_sub,dof_int,dof_int_L,dof_int_R
   integer :: i,j,k,l,ni,nb,ndof,noff
   integer :: isub,jsub,kib,kbb,nnz_sub,nnz_ib
   VTYPE   :: x
!
!..sparse MKL
   integer :: mkl_stat
   integer, allocatable :: Aib_row(:), Aib_col(:)
   VTYPE  , allocatable :: Aib_val(:)
   type (SPARSE_MATRIX_T) :: Aib_sparse
   type (MATRIX_DESCR)    :: Aib_descr
   VTYPE, allocatable :: Aii(:,:),Bi(:)
!
!..timer
   real(8) :: time_stamp,start_time
!
!..info (verbose output if true)
   logical :: info = .false.
   logical :: dbg  = .false. ! print additional info for debugging purposes
!
!..auxiliary storage for reallocation
   integer, allocatable :: tmp_irn(:),tmp_jcn(:)
   VTYPE  , allocatable :: tmp_val(:)
!
! -----------------------------------------------------------------------
!
 123 format('[',I3,'] ',A10,': ',I8)
!
   call MPI_COMM_RANK(mumps%COMM, mRANK ,ierr)
   call MPI_COMM_SIZE(mumps%COMM, mPROCS,ierr)
   if (mPROCS .ne. nproc) then
      write(*,*) 'par_fiber: mPROCS != nproc : ',mPROCS,nproc
   endif
   if (dbg) write(*,123) RANK,'mRANK',mRANK
   if (dbg) write(*,123) RANK,'mPROCS',mPROCS
!
!..if the current problem is "small enough", solve with parallel MUMPS
   if (mPROCS .le. mSUB_PROCS) then
      call par_solve(mumps)
      goto 90
   endif
!
   if (mRANK .eq. ROOT) then
      write(*,5010) '[',RANK,'] par_fiber (level ',level,'): mPROCS = ', mPROCS, &
                                                      ', mSUB_PROCS = ', mSUB_PROCS
      !write(*,*) ' - solving distributed sparse problem: mPROCS     = ',mPROCS
      !write(*,*) ' - splitting into subproblems of size: mSUB_PROCS = ',mSUB_PROCS
 5010 format(A,I4,A,I2,A,I4,A,I4)
   endif
!
!..currently, we assume that we are doing nested dissection on 2^k processors
   if (MOD(mPROCS,mSUB_PROCS) .ne. 0) then
      if (mRANK .eq. ROOT) write(*,*) 'par_fiber: mPROCS != 2^k. Solving without nested dissection.'
      call par_solve(mumps)
      goto 90
   endif
!
   call MPI_BARRIER(mumps%COMM, ierr); start_time = MPI_Wtime()
!
!..divide into "small enough" subproblems (groups) and solve separately
!  [Step 1: define communicators for the subproblems]
   group = mRANK / mSUB_PROCS
   color = group
   key   = MOD(mRANK,mSUB_PROCS)
   call MPI_COMM_SPLIT(mumps%COMM,color,key, mpi_comm_sub,ierr)
   call MPI_COMM_RANK(mpi_comm_sub, mSUB_RANK,ierr)
   if (mSUB_RANK .ne. key) then
      write(*,5050) 'RANK, mSUB_RANK, key = ',RANK,mSUB_RANK,key
 5050 format(A,I4,',',I4,',',I4)
   endif
   call mumps_start(mumps_sub,mpi_comm_sub)
!
!..and the interface problem (needs one host from each group)
!  [Step 2: define communicator for the separator problem]
   color = MPI_UNDEFINED; mINT_RANK = -1
   key   = group
   if (mSUB_RANK .eq. ROOT) color = 0
   call MPI_COMM_SPLIT(mumps%COMM,color,key, mpi_comm_int,ierr)
   if (mSUB_RANK .eq. ROOT) then
      call MPI_COMM_RANK(mpi_comm_int, mINT_RANK ,ierr)
      call MPI_COMM_SIZE(mpi_comm_int, mINT_PROCS,ierr)
      if (mINT_RANK .ne. key) then
         write(*,5050) 'RANK, mINT_RANK, key = ',RANK,mINT_RANK,key
      endif
      call mumps_start(mumps_int,mpi_comm_int)
   endif
   if (dbg) write(*,123) RANK,'mSUB_RANK',mSUB_RANK
   if (dbg) write(*,123) RANK,'mINT_RANK',mINT_RANK
   if (dbg) write(*,123) RANK,'nrdof',nrdof(mRANK+1)
!
!..calculate number of (interior) dofs for the subproblem
   dof_sub = sum(nrdof(group*mSUB_PROCS+1:group*mSUB_PROCS+mSUB_PROCS-1))
   if (dbg) write(*,123) RANK,'dof_sub',dof_sub
!
!..calculate number of (interface) dofs for the subproblem (both left and right interface)
   dof_int_L = 0
   if (group .gt. 0) dof_int_L = nrdof(group*mSUB_PROCS)
   dof_int_R = nrdof(group*mSUB_PROCS + mSUB_PROCS)
   dof_int = dof_int_L + dof_int_R
   if (dbg) write(*,123) RANK,'dof_int',dof_int
   if (dbg) write(*,123) RANK,'dof_int_L',dof_int_L
   if (dbg) write(*,123) RANK,'dof_int_R',dof_int_R
!
!..calculate dof offset for a group
   dof_off = sum(nrdof(1:group*mSUB_PROCS))
   if (dbg) write(*,123) RANK,'dof_off',dof_off
!
!  in the input MUMPS instance, the global RHS is only valid on (mRank.eq.ROOT) processor
!  even though the other procs do have memory allocated
!  --> distribute RHS to group reps (mSUB_RANK.eq.0) who can then build
!      the subproblem RHS from it (can be optimized later, maybe via scatter...)
!
!..broadcast global RHS from host (group rep of group 0, i.e., mINT_RANK=0)
!  to group reps of other groups
   if (mSUB_RANK .eq. ROOT) then
      count = mumps%N; src = ROOT
      call MPI_BCAST(mumps%RHS,count,MPI_VTYPE,src,mumps_int%COMM,ierr)
   endif
!
!..allocate subproblem data structures
   mumps_sub%N = dof_sub
   nnz_sub = int(mumps%NNZ_loc,4) ! -> TODO: can we give a better estimate (upper bound)
   mumps_sub%NNZ_loc = nnz_sub    ! -> adjust later before calling solver
   mumps_sub%NRHS = NR_RHS + dof_int
   mumps_sub%LRHS = dof_sub
   allocate(mumps_sub%IRN_loc(nnz_sub))
   allocate(mumps_sub%JCN_loc(nnz_sub))
   allocate(mumps_sub%A_loc(nnz_sub))
   allocate(mumps_sub%RHS(dof_sub * (NR_RHS+dof_int)))
   mumps_sub%RHS = ZERO
!
!..interface contributions should only come from the first and last proc
   if ((mSUB_RANK.eq.ROOT) .or. (mSUB_RANK.eq.mSUB_PROCS-1)) then
!  ...create dense interface stiffness matrix Aii and load Bi
      allocate(Aii(dof_int,dof_int)); Aii = ZERO
      if (mSUB_RANK.eq.ROOT) then
         allocate(Bi(dof_int)); Bi = ZERO
      endif
!
!  ...create sparse matrix Aib in COO format
      nnz_ib = int(mumps%NNZ_loc,4) ! -> TODO: can we give a better estimate (upper bound)
      allocate(Aib_val(nnz_ib))
      allocate(Aib_row(nnz_ib))
      allocate(Aib_col(nnz_ib))
   endif
!
!..create matrices and rhs for the subproblems
!  [Steps 3-4: Assembling the subproblems and the separator problem]
   if (mSUB_RANK .eq. ROOT) then
!  ...split global RHS into interface RHS vector and subproblem RHS vector
      do j=0,NR_RHS-1
!     ...1. assemble RHS vector for interface (separator) problem
!           each group rep does this only for their left side (lower indices)
!           to avoid duplicating already assembled RHS interface contributions
!           coming from the group rep on the right side (higher indices)
         do i=1,dof_int_L
            ! add to Bi
            Bi(j*dof_int + i) = mumps%RHS(j*mumps%N + dof_off-dof_int_L+i)
         enddo
!     ...2. assemble RHS vector for subproblem
         do i=1,dof_sub
            ! add to Bsub
            mumps_sub%RHS(j*dof_sub + i) = mumps%RHS(j*mumps%N + dof_off+i)
         enddo
      enddo
   endif
   kbb = 0; kib = 0 ! nnz counters
   do k = 1,int(mumps%NNZ_loc,4)
      i = mumps%IRN_loc(k)
      j = mumps%JCN_loc(k)
      x = mumps%A_loc(k)
      if (i - dof_off .le. 0) then ! i = interface dof (row) to the left side (verify rank)
         if (mSUB_RANK .ne. ROOT) then ! CHECK
            write(*,*) 'par_fiber: mSUB_RANK .ne. ROOT',mSUB_RANK,ROOT; stop
         endif
         if (i-dof_off+dof_int_L .le. 0) then ! CHECK
            write(*,*) 'par_fiber: unexpected coupling. stop.'; stop
         endif
         if (j - dof_off .le. 0) then ! j = interface dof (col) to the left side
            if (mSUB_RANK .ne. ROOT) then ! CHECK
               write(*,*) 'par_fiber: mSUB_RANK .ne. ROOT',mSUB_RANK,ROOT; stop
            endif
            if (j-dof_off+dof_int_L .le. 0) then ! CHECK
               write(*,*) 'par_fiber: unexpected coupling. stop.'; stop
            endif
!        ...add to A_ii
            isub = i - dof_off + dof_int_L
            jsub = j - dof_off + dof_int_L
            Aii(isub,jsub) = Aii(isub,jsub) + x
         elseif (j-dof_off .gt. dof_sub) then ! CHECK: left most and right most should not be connected
            write(*,*) 'par_fiber: unexpected coupling. stop.'; stop
         else ! j = bubble dof (col)
!        ...add to A_ib
            kib = kib + 1
            if (kib > nnz_ib) then ! CHECK
               write(*,*) 'kib > nnz_ib',kib,nnz_ib; stop
            endif
            isub = i - dof_off + dof_int_L
            jsub = j - dof_off
            Aib_row(kib) = isub
            Aib_col(kib) = jsub
            Aib_val(kib) = x
         endif
!
      elseif (i - dof_off .gt. dof_sub) then ! interface dof to the right side (verify rank)
         if (mSUB_RANK .ne. mSUB_PROCS-1) then
            write(*,*) 'par_fiber: mSUB_RANK .ne. mSUB_PROCS-1',mSUB_RANK,mSUB_PROCS-1; stop
         endif
         if (i-dof_off .gt. dof_sub+dof_int_R) then ! CHECK
            write(*,*) 'par_fiber: unexpected coupling. stop.'; stop
         endif
         if (j-dof_off .le. 0) then ! CHECK: left most and right most should not be connected
            write(*,*) 'par_fiber: unexpected coupling. stop.'; stop
         elseif (j - dof_off .gt. dof_sub) then ! CHECK
            if (mSUB_RANK .ne. mSUB_PROCS-1) then
               write(*,*) 'par_fiber: mSUB_RANK .ne. mSUB_PROCS-1',mSUB_RANK,mSUB_PROCS-1; stop
            endif
            if (j-dof_off .gt. dof_sub+dof_int_R) then ! CHECK
               write(*,*) 'par_fiber: unexpected coupling. stop.'; stop
            endif
!        ...add to A_ii
            isub = i - (dof_off+dof_sub) + dof_int_L
            jsub = j - (dof_off+dof_sub) + dof_int_L
            Aii(isub,jsub) = Aii(isub,jsub) + x
         else ! j = bubble dof (col)
!        ...add to A_ib
            kib = kib + 1
            if (kib > nnz_ib) then ! CHECK
               write(*,*) 'kib > nnz_ib',kib,nnz_ib; stop
            endif
            isub = i - (dof_off+dof_sub) + dof_int_L
            jsub = j - dof_off
            Aib_row(kib) = isub
            Aib_col(kib) = jsub
            Aib_val(kib) = x
         endif
!
!  ...i = bubble dof (row)
      else
!        [every proc except left most and right most in group should only enter this part of the branch]
         if (j - dof_off .le. 0) then ! j = interface dof (col) to the left side
            if (mSUB_RANK .ne. ROOT) then ! CHECK
               write(*,*) 'par_fiber: mSUB_RANK .ne. ROOT',mSUB_RANK,ROOT; stop
            endif
            if (j-dof_off+dof_int_L .le. 0) then ! CHECK
               write(*,*) 'par_fiber: unexpected coupling. stop.'; stop
            endif
!        ...add to A_bi
            isub = i - dof_off
            jsub = j - dof_off + dof_int_L
            l = NR_RHS*dof_sub + (jsub-1)*dof_sub + isub
            mumps_sub%RHS(l) = mumps_sub%RHS(l) + x
         elseif (j - dof_off .gt. dof_sub) then ! j = interface dof (col) to the right side
            if (mSUB_RANK .ne. mSUB_PROCS-1) then ! CHECK
               write(*,*) 'par_fiber: mSUB_RANK .ne. mSUB_PROCS-1',mSUB_RANK,mSUB_PROCS-1; stop
            endif
            if (j-dof_off .gt. dof_sub+dof_int_R) then ! CHECK
               write(*,*) 'par_fiber: unexpected coupling. stop.'; stop
            endif
!        ...add to A_bi
            isub = i - dof_off
            jsub = j - (dof_off+dof_sub) + dof_int_L
            l = NR_RHS*dof_sub + (jsub-1)*dof_sub + isub
            mumps_sub%RHS(l) = mumps_sub%RHS(l) + x
         else ! j = bubble dof (col)
!           [every proc except left most and right most should only enter this part of the branch]
!        ...add to A_bb
            kbb = kbb + 1
            if (kbb > nnz_sub) then ! CHECK
               write(*,*) 'kbb > nnz_sub',kbb,nnz_sub; stop
            endif
            isub = i - dof_off
            jsub = j - dof_off
            mumps_sub%IRN_loc(kbb) = isub
            mumps_sub%JCN_loc(kbb) = jsub
            mumps_sub%A_loc(kbb) = x
         endif
      endif
   enddo
!
   deallocate(mumps%A_loc,mumps%IRN_loc,mumps%JCN_loc)
!
   if (dbg) write(*,123) RANK,'kib',kib
   if (dbg) write(*,123) RANK,'kbb',kbb
!
!..verify indices are valid
   if ((mSUB_RANK.eq.ROOT) .or. (mSUB_RANK.eq.mSUB_PROCS-1)) then
      do k=1,kib
         i = Aib_row(k)
         j = Aib_col(k)
         if ((i .le. 0) .or. (i .gt. dof_int)) then
            write(*,*) 'par_fiber: invalid index: i = ', i; stop
         endif
         if ((j .le. 0) .or. (j .gt. dof_sub)) then
            write(*,*) 'par_fiber: invalid index: j = ', j; stop
         endif
      enddo
   endif
   do k=1,kbb
      i = mumps_sub%IRN_loc(k)
      j = mumps_sub%JCN_loc(k)
      if ((i .le. 0) .or. (i .gt. dof_sub)) then
         write(*,*) 'par_fiber: invalid index: i = ', i; stop
      endif
      if ((j .le. 0) .or. (j .gt. dof_sub)) then
         write(*,*) 'par_fiber: invalid index: j = ', j; stop
      endif
   enddo
!
!..set non-zero counters
   nnz_sub = kbb
   mumps_sub%NNZ_loc = nnz_sub
!
!  TODO: Check if reallocating the arrays to the right size changes anything
!   allocate(tmp_irn(nnz_sub),tmp_jcn(nnz_sub),tmp_val(nnz_sub))
!   tmp_irn(1:nnz_sub) = mumps_sub%IRN_loc(1:nnz_sub)
!   tmp_jcn(1:nnz_sub) = mumps_sub%JCN_loc(1:nnz_sub)
!   tmp_val(1:nnz_sub) = mumps_sub%A_loc  (1:nnz_sub)
!   deallocate(mumps_sub%IRN_loc,mumps_sub%JCN_loc,mumps_sub%A_loc)
!   allocate(mumps_sub%IRN_loc(nnz_sub),mumps_sub%JCN_loc(nnz_sub),mumps_sub%A_loc(nnz_sub))
!   mumps_sub%IRN_loc(1:nnz_sub) = tmp_irn(1:nnz_sub); deallocate(tmp_irn)
!   mumps_sub%JCN_loc(1:nnz_sub) = tmp_jcn(1:nnz_sub); deallocate(tmp_jcn)
!   mumps_sub%A_loc  (1:nnz_sub) = tmp_val(1:nnz_sub); deallocate(tmp_val)
!
!  gather (reduce) RHS information on host [b_sub | A_bi]
!  b_sub is already on the host since it was gathered from global RHS in mumps instance
!  ...but we still need to collect A_bi information (this should only come from last proc)
!  ...the corresponding A_bi values are in columns referring to interface dof on the right side,
!     hence we can overwrite the values (which should currently be zeros on the host)
   if (mSUB_RANK .eq. ROOT) then
      src = mSUB_PROCS-1; count = dof_sub*dof_int_R; tag = 1
      call MPI_RECV(mumps_sub%RHS(dof_sub*(NR_RHS+dof_int_L)+1:dof_sub*(NR_RHS+dof_int)),count,MPI_VTYPE,src,tag,mumps_sub%COMM, MPI_STATUS_IGNORE,ierr)
   elseif (mSUB_RANK .eq. mSUB_PROCS-1) then
      rcv = ROOT; count = dof_sub*dof_int_R; tag = 1
      call MPI_SEND(mumps_sub%RHS(dof_sub*(NR_RHS+dof_int_L)+1:dof_sub*(NR_RHS+dof_int)),count,MPI_VTYPE,rcv,tag,mumps_sub%COMM, ierr)
   endif
!
!..percentage increase in estimated workspace for global interface problem
   mumps_sub%icntl(14) = 30
!
!..solve subproblems
   call par_solve(mumps_sub)
!
   if (associated(mumps_sub%A_loc))   deallocate(mumps_sub%A_loc)
   if (associated(mumps_sub%IRN_loc)) deallocate(mumps_sub%IRN_loc)
   if (associated(mumps_sub%JCN_loc)) deallocate(mumps_sub%JCN_loc)
!
!..if proc is not participating in solving the interface problem, then skip the following section
!  the right proc of each group skips later, it must still communicate interactions to the host
   if ((mSUB_RANK.ne.ROOT) .and. (mSUB_RANK.ne.mSUB_PROCS-1)) goto 60
!
!..send A_ib matrix to host (left proc) to be assembled
   if (mSUB_RANK .eq. ROOT) then
      src = mSUB_PROCS-1; count = 1; tag = 1
      call MPI_RECV(l,count,MPI_INTEGER,src,tag,mumps_sub%COMM, MPI_STATUS_IGNORE,ierr)
      if (kib + l > nnz_ib) then ! CHECK (initial allocation was insufficient)
         write(*,*) 'MUST REALLOCATE since kib + l > nnz_ib',kib,l,nnz_ib
         allocate(tmp_irn(kib + l),tmp_jcn(kib + l),tmp_val(kib + l))
         tmp_irn(1:kib) = Aib_row(1:kib); call move_alloc(tmp_irn, Aib_row)
         tmp_jcn(1:kib) = Aib_col(1:kib); call move_alloc(tmp_jcn, Aib_col)
         tmp_val(1:kib) = Aib_val(1:kib); call move_alloc(tmp_val, Aib_val)
      endif
      count = l; tag = 2
      call MPI_RECV(Aib_row(kib+1:kib+l),count,MPI_INTEGER,src,tag,mumps_sub%COMM, MPI_STATUS_IGNORE,ierr)
      tag = 3
      call MPI_RECV(Aib_col(kib+1:kib+l),count,MPI_INTEGER,src,tag,mumps_sub%COMM, MPI_STATUS_IGNORE,ierr)
      tag = 4
      call MPI_RECV(Aib_val(kib+1:kib+l),count,MPI_VTYPE,src,tag,mumps_sub%COMM, MPI_STATUS_IGNORE,ierr)
      nnz_ib = kib + l
      ! reduction (sum) over Aii (dense matrix)
      count = dof_int*dof_int_R; tag = 5
      call MPI_RECV(Aii(1:dof_int,dof_int_L+1:dof_int),count,MPI_VTYPE,src,tag,mumps_sub%COMM, MPI_STATUS_IGNORE,ierr)
   elseif (mSUB_RANK .eq. mSUB_PROCS-1) then
      rcv = ROOT; count = 1; tag = 1
      call MPI_SEND(kib,count,MPI_INTEGER,rcv,tag,mumps_sub%COMM, ierr)
      count = kib; tag = 2
      call MPI_SEND(Aib_row(1:kib),count,MPI_INTEGER,rcv,tag,mumps_sub%COMM, ierr)
      tag = 3
      call MPI_SEND(Aib_col(1:kib),count,MPI_INTEGER,rcv,tag,mumps_sub%COMM, ierr)
      tag = 4
      call MPI_SEND(Aib_val(1:kib),count,MPI_VTYPE,rcv,tag,mumps_sub%COMM, ierr)
      ! reduction (sum) over Aii (dense matrix)
      count = dof_int*dof_int_R; tag = 5
      call MPI_SEND(Aii(1:dof_int,dof_int_L+1:dof_int),count,MPI_VTYPE,rcv,tag,mumps_sub%COMM, ierr)
      deallocate(Aib_val,Aib_row,Aib_col,Aii)
      goto 60
   endif
!
!..create matrix and rhs for the interface problem
!  each group (representative) works on assembling their part of the interface problem
   ni = dof_int; nb = dof_sub
!
   call MPI_BARRIER(mumps_int%COMM, ierr); time_stamp = MPI_Wtime()
!
   Aib_descr%type = SPARSE_MATRIX_TYPE_GENERAL
#if C_MODE
   mkl_stat = MKL_SPARSE_Z_CREATE_COO(Aib_sparse,              &
                                      SPARSE_INDEX_BASE_ONE,   &
                                      ni,nb,nnz_ib,            &
                                      Aib_row(1:nnz_ib),Aib_col(1:nnz_ib),Aib_val(1:nnz_ib))
   mkl_stat = MKL_SPARSE_Z_MM(SPARSE_OPERATION_NON_TRANSPOSE,              &
                              -ZONE,                                       &
                              Aib_sparse,Aib_descr,                        &
                              SPARSE_LAYOUT_COLUMN_MAJOR,                  &
                              mumps_sub%RHS(1:nb*NR_RHS),                  &
                              NR_RHS, nb, ZONE, Bi, ni)
   mkl_stat = MKL_SPARSE_Z_MM(SPARSE_OPERATION_NON_TRANSPOSE,              &
                              -ZONE,                                       &
                              Aib_sparse,Aib_descr,                        &
                              SPARSE_LAYOUT_COLUMN_MAJOR,                  &
                              mumps_sub%RHS(nb*NR_RHS+1:nb*(NR_RHS+ni)),   &
                              ni, nb, ZONE, Aii, ni)
#else
   mkl_stat = MKL_SPARSE_D_CREATE_COO(Aib_sparse,              &
                                      SPARSE_INDEX_BASE_ONE,   &
                                      ni,nb,nnz_ib,            &
                                      Aib_row(1:nnz_ib),Aib_col(1:nnz_ib),Aib_val(1:nnz_ib))
   mkl_stat = MKL_SPARSE_D_MM(SPARSE_OPERATION_NON_TRANSPOSE,              &
                              -ZONE,                                       &
                              Aib_sparse,Aib_descr,                        &
                              SPARSE_LAYOUT_COLUMN_MAJOR,                  &
                              mumps_sub%RHS(1:nb*NR_RHS),                  &
                              NR_RHS, nb, ZONE, Bi, ni)
   mkl_stat = MKL_SPARSE_D_MM(SPARSE_OPERATION_NON_TRANSPOSE,              &
                              -ZONE,                                       &
                              Aib_sparse,Aib_descr,                        &
                              SPARSE_LAYOUT_COLUMN_MAJOR,                  &
                              mumps_sub%RHS(nb*NR_RHS+1:nb*(NR_RHS+ni)),   &
                              ni, nb, ZONE, Aii, ni)
#endif
!
   call MPI_BARRIER(mumps_int%COMM, ierr)
   if (mINT_RANK.eq.ROOT .and. info) then
      write(*,7891) '[',RANK,'] INTERFACE ASSEMBLY Sparse MKL: ', MPI_Wtime()-time_stamp
7891  format(A,I4,A,F12.5)
   endif
!
   mkl_stat = MKL_SPARSE_DESTROY(Aib_sparse)
   deallocate(Aib_val,Aib_row,Aib_col)
!
!..percentage increase in estimated workspace for global interface problem
   mumps_int%icntl(14) = 30
!
!..memory allocation for global MUMPS solve (remaining coupled interface problem)
   ndof = sum(nrdof(mSUB_PROCS:mPROCS : mSUB_PROCS)) ! total number of interface dof
   if (dbg) write(*,123) RANK,'ndof',ndof
   mumps_int%N = ndof ! total number of interface dof
   mumps_int%NNZ_loc = dof_int * dof_int ! the local matrix is in fact densely populated
   mumps_int%NRHS = NR_RHS ! currently assumed one right-hand side only
   ! mumps_int%LRHS = ndof
   allocate(mumps_int%IRN_loc(mumps_int%NNZ_loc))
   allocate(mumps_int%JCN_loc(mumps_int%NNZ_loc))
   allocate(mumps_int%A_loc(mumps_int%NNZ_loc))
   allocate(mumps_int%RHS(mumps_int%N * NR_RHS))
   mumps_int%RHS=ZERO
!
!..convert dense local interface matrix into sparse global interface matrix
!  and convert dense local load vector into dense global load vector; currently assuming NR_RHS=1
!  calculate offset of this group's interface dof in the global interface solution
   noff = sum(nrdof(mSUB_PROCS : mRANK-1 : mSUB_PROCS))
   if (dbg) write(*,123) RANK,'noff',noff
!
!$OMP PARALLEL DO PRIVATE(i,j,l,isub,jsub)
   do jsub = 1,dof_int ! col
      j = noff+jsub
      mumps_int%RHS(j) = Bi(jsub)
      do isub = 1,dof_int ! row
         i = noff+isub
         l = (jsub-1)*dof_int + isub
         mumps_int%A_loc(l) = Aii(isub,jsub)
         mumps_int%IRN_loc(l) = i
         mumps_int%JCN_loc(l) = j
      enddo
   enddo
!$OMP END PARALLEL DO
!
   deallocate(Aii,Bi)
!
!..gather RHS vector information on host
   count = mumps_int%N
   if (RANK .eq. ROOT) then
      call MPI_REDUCE(MPI_IN_PLACE,mumps_int%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_int%COMM, ierr)
   else
      call MPI_REDUCE(mumps_int%RHS,mumps_int%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_int%COMM, ierr)
   endif
!
!..solve interface problem
   if (RECURSIVE_NESTED) then
      nrdof_int(1:mINT_PROCS) = nrdof(mSUB_PROCS:mPROCS : mSUB_PROCS)
      call par_fiber(mumps_int,nrdof_int,mINT_PROCS,level+1)
   else
      call par_solve(mumps_int)
   endif
!
   if (associated(mumps_int%A_loc))   deallocate(mumps_int%A_loc)
   if (associated(mumps_int%IRN_loc)) deallocate(mumps_int%IRN_loc)
   if (associated(mumps_int%JCN_loc)) deallocate(mumps_int%JCN_loc)
!
!..backsubstitute to obtain subproblem solutions
!  1. broadcast interface solution to group representatives
   count = mumps_int%N; src = ROOT
   call MPI_BCAST(mumps_int%RHS,count,MPI_VTYPE,src,mumps_int%COMM, ierr)
!
!..2. Perform backsubstitution to obtain subproblem bubble solution (assuming NR_RHS=1 currently)
!     NOTE: All vectors and matrices here are dense
!     Compute (dense)  vector: x_b =  -[A_bb^-1 A_bi] x_i + [A_bb^-1 * l_b]
!     store solution in mumps_sub%RHS(1:mumps_sub%N*NR_RHS)
   ni = dof_int; nb = dof_sub
   ! calculate offset of this group's interface dof in the global interface solution
   noff = sum(nrdof(mSUB_PROCS : mRANK-1 : mSUB_PROCS))
   call stc_bwd(ni,nb,                                      &
                mumps_sub%RHS(nb*NR_RHS+1:nb*(NR_RHS+ni)),  &
                mumps_int%RHS(noff+1:noff+ni), mumps_sub%RHS(1:nb*NR_RHS))
!
!..assemble back into global solution vector (on host); currently assuming NR_RHS=1
!  mumps%RHS should hold the global solution in the end on the host proc (also the host on mumps_int)
!  every group rep writes the subproblem solution into mumps%RHS for the dofs that the group owns
!  then we reduce over mumps%RHS to collect all values
   mumps%RHS = ZERO
!  insert subproblem solution into mumps%RHS for group owned dofs
   mumps%RHS(dof_off+1        :dof_off+dof_sub          ) = mumps_sub%RHS(1:dof_sub)
   mumps%RHS(dof_off+dof_sub+1:dof_off+dof_sub+dof_int_R) = mumps_int%RHS(noff+dof_int_L+1:noff+dof_int)
!
   count = mumps%N
   if (mINT_RANK .eq. ROOT) then
      call MPI_REDUCE(MPI_IN_PLACE,mumps%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_int%COMM, ierr)
   else
      call MPI_REDUCE(mumps%RHS,mumps%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_int%COMM, ierr)
   endif
!
!..cleanup
   call mumps_destroy(mumps_int)
!
 60 continue
!
   call MPI_BARRIER(mumps%COMM, ierr)
   if (mRANK.eq.ROOT) then
      time_stamp = MPI_Wtime()-start_time
      write(*,9090) '[',RANK,'] par_fiber (level ',level,'): ',time_stamp,' seconds'
 9090 format(A,I4,A,I2,A,F12.5,A)
   endif
!
!..cleanup
   call mumps_destroy(mumps_sub)
!
 90 continue
!
end subroutine par_fiber

#else

subroutine par_fiber(mumps,nrdof,nproc,level)
   use mpi_param , only: RANK,ROOT
   use par_mumps
   implicit none
#if C_MODE
   type (ZMUMPS_STRUC), intent(inout) :: mumps
#else
   type (DMUMPS_STRUC), intent(inout) :: mumps
#endif
   integer, intent(in) :: nproc
   integer, intent(in) :: nrdof(nproc)
   integer, intent(in) :: level
   if (RANK .eq. ROOT) then
      write(*,*) 'par_fiber: HP3D_USE_INTEL_MKL = 0. Dependency is required.'
   endif
end subroutine par_fiber

#endif
