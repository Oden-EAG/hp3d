!
#include "implicit_none.h"
! -----------------------------------------------------------------------
!
!    routine name       - par_mumps_sc
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 2019
!
!    purpose            - interface for distributed MUMPS solver
!                       - routine computes global stiffness matrix
!                         and global load vector and solves the global
!                         linear system with MUMPS
!                       - the assembly is computed in parallel with one
!                         MPI process per subdomain, and OpenMP
!                         parallelization within each subdomain
!                       - this routine supports both computation with or
!                         without static condensation (uses module stc)
!
!    in                 - mtype: 'H':  Hermitian (for complex)/
!                                      Symmetric (for real)
!                                      case for Lapack routines
!                                'G':  General case for Lapack routines
!               note: MUMPS does currently not support LU factorization
!                     of Hermitian matrices. In that case, it will use
!                     the non-symmetric matrix setting
!
! -----------------------------------------------------------------------
subroutine par_mumps_sc(mtype)
!
   use data_structure3D, only: NRNODS, NRELES, get_subd
   use assembly,         only: NR_RHS, MAXDOFM, MAXDOFS,       &
                               MAXbrickH, MAXmdlbH, NRHVAR,    &
                               MAXbrickE, MAXmdlbE, NREVAR,    &
                               MAXbrickV, MAXmdlbV, NRVVAR,    &
                               MAXbrickQ, NRQVAR,              &
                               NEXTRACT, IDBC, ZDOFD, ZERO,    &
                               ALOC, BLOC, AAUX, ZAMOD, ZBMOD, &
                               NR_PHYSA, MAXNODM
   use assembly_sc
   use control,   only: ISTC_FLAG
   use stc,       only: HERM_STC,CLOC,stc_alloc,stc_dealloc,stc_get_nrdof
   use par_mumps, only: mumps_par,mumps_start,mumps_destroy
   use par_mesh , only: DISTRIBUTED,HOST_MESH
   use mpi_param, only: RANK,ROOT
   use MPI      , only: MPI_SUM,MPI_REAL8,MPI_COMPLEX16
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
   integer    :: nrdof_H,nrdof_E,nrdof_V,nrdof_Q
   integer    :: nrdofm,nrdofc,nrnodm,nrdof,nrdof_mdl,ndof
   integer    :: iel,mdle,subd,i,j,k,l,k1,k2,nod,idec
   integer(8) :: nnz,nnz_loc
!
!..MPI variables
   integer :: count,src,ierr
!
!..dummy variables
   integer :: nvoid
   VTYPE   :: zvoid
! 
!..workspace for celem
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
!
   integer    :: mdle_list(NRELES)
   integer(8) :: elem_nnz_loc(NRELES)
!
!..workspace for right-hand side and solution vector
   VTYPE, allocatable :: RHS(:),SOL(:)
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time,time_stamp
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
   if ((.not. DISTRIBUTED) .or. HOST_MESH) then
      if (RANK .eq. ROOT) then
         write(*,*) 'par_mumps_sc: mesh is not distributed (or on host).'
         write(*,*) 'calling mumps_sc from host...'
         call mumps_sc(mtype)
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
1000  format(' par_mumps_sc: STARTED')
      write(*,*)
   endif
!
!..TODO multiple right-hand sides
   NR_RHS = 1
   call mumps_start
!
!..calculate maximum number of dofs for a modified element
!  with static condensation
   if (ISTC_FLAG) then
      MAXDOFM = (MAXbrickH-MAXmdlbH)*NRHVAR   &
              + (MAXbrickE-MAXmdlbE)*NREVAR   &
              + (MAXbrickV-MAXmdlbV)*NRVVAR
!  without static condensation
   else
      MAXDOFM = MAXbrickH*NRHVAR    & 
              + MAXbrickE*NREVAR    &
              + MAXbrickV*NRVVAR    &
              + MAXbrickQ*NRQVAR
   endif
!
! ----------------------------------------------------------------------
!  STEP 1 : 1ST LOOP THROUGH ELEMENTS, 1ST CALL TO CELEM TO GET INFO
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      if (RANK .eq. ROOT) write(*,1001)
 1001 format(' STEP 1 started : Get assembly info')
      call MPI_BARRIER(mumps_par%COMM, ierr)
      start_time = MPI_Wtime()
   endif
!
   allocate(MAXDOFS(NR_PHYSA))
   MAXDOFS = 0; MAXDOFM = 0
!
!..allocate and initialize offsets
   allocate(NFIRSTH(NRNODS)); NFIRSTH = -1
   allocate(NFIRSTE(NRNODS)); NFIRSTE = -1
   allocate(NFIRSTV(NRNODS)); NFIRSTV = -1
   if (.not. ISTC_FLAG) then
      allocate(NFIRSTQ(NRNODS)); NFIRSTQ = -1
   endif
!
!..Step 1: determine the first dof offsets for active nodes
   nrdof_H = 0; nrdof_E = 0; nrdof_V = 0; nrdof_Q = 0; nrdof_mdl = 0
!
   mdle = 0; idec = 1
!..non-zero counters for element offsets in distributed sparse stiffness matrix
!  using 64 bit integers nnz and nnz_loc
   nnz     = 0_8; ! global number of non-zeros in matrix (counts duplicate indices)
   nnz_loc = 0_8; ! local  number of non-zeros in matrix (counts duplicate indices)
   elem_nnz_loc(1:NRELES) = 0_8 ! local element offsets for subdomain matrix
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(iel) = mdle
      call get_subd(mdle, subd)
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
!  ...k      = number of non-zeros entries in element stiffness matrix
      k = nrdofc**2
!
!  ...global counters for OpenMP
      nnz = nnz + int8(k)
!
      if (subd .eq. RANK) then
!     ...subdomain counters for OpenMP
         elem_nnz_loc(iel) = nnz_loc
         nnz_loc = nnz_loc + int8(k)
!
!     ...update the maximum number of local dof
         do i=1,NR_PHYSA
            MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
         enddo
!     ...update the maximum number of modified element dof in the expanded mode
         MAXDOFM = max0(MAXDOFM,nrdofm)
      endif
!
!  ...compute offsets for H1 dof
      do i = 1,nrnodm
         nod = nodm(i)
!     ...avoid repetition
         if (NFIRSTH(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTH(nod) = nrdof_H
!     ...update the H1 dof counter
         nrdof_H = nrdof_H + ndofmH(i)
      enddo
!
!  ...compute offsets for H(curl) dof
      do i = 1,nrnodm
         nod = nodm(i)
!     ...avoid repetition
         if (NFIRSTE(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTE(nod) = nrdof_E
!     ...update the H(curl) dof counter
         nrdof_E = nrdof_E + ndofmE(i)
      enddo
!
!  ...compute offsets for H(div) dof
      do i=1,nrnodm
         nod = nodm(i)
!     ...avoid repetition
         if (NFIRSTV(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTV(nod) = nrdof_V
!     ...update the H(div) dof counter
         nrdof_V = nrdof_V + ndofmV(i)
      enddo
!
      if (.not. ISTC_FLAG) then
         nod = nodm(nrnodm)
!     ...avoid repetition
         if (NFIRSTQ(nod).lt.0) then
!        ...store the first dof offset
            NFIRSTQ(nod) = nrdof_Q
!        ...update the L2 dof counter
            nrdof_Q = nrdof_Q + ndofmQ(nrnodm)
         endif
      endif
!
!  ...compute number of bubble dof (nrdof_mdl)
      if (ISTC_FLAG) then
         call stc_get_nrdof(mdle, nrdofi,nrdofb)
         nrdof_mdl = nrdof_mdl + sum(nrdofb)
      endif
!
!..end of loop through elements
   enddo
!
!..total number of (interface) dof is nrdof
   nrdof = nrdof_H + nrdof_E + nrdof_V + nrdof_Q
!
   NRDOF_CON = nrdof
   NRDOF_TOT = nrdof + nrdof_mdl
!
   if (nrdof .eq. 0) then
      deallocate(MAXDOFS,NFIRSTH,NFIRSTE,NFIRSTV)
      if (.not. ISTC_FLAG) deallocate(NFIRSTQ)
      if (RANK .eq. ROOT) write(*,*) 'par_mumps_sc: nrdof = 0. returning.'
      return
   endif
!
! ----------------------------------------------------------------------
!  END OF STEP 1
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
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
   call MPI_BARRIER(mumps_par%COMM, ierr)
   if (RANK .eq. ROOT) then
      write(*,2010) '[', RANK, '] Number of dof  : nrdof_con = ', NRDOF_CON
      write(*,2010) '[', RANK, ']                  nrdof_tot = ', NRDOF_TOT
      write(*,2010) '[', RANK, '] Total non-zeros: nnz       = ', nnz
   endif
   write(*,2010) '[', RANK, '] Local non-zeros: nnz_loc   = ', nnz_loc
2010 format(A,I4,A,I12)
!
!..use 64bit parallel analysis if nnz_loc > 2B
!  (sequential metis/scotch using 32bit currently)
   if (nnz > 2e9) mumps_par%icntl(28) = 2
!
!..memory allocation for load assembly
   allocate(RHS(nrdof)); RHS=ZERO
!
!..memory allocation for MUMPS solver
   mumps_par%N = nrdof
   mumps_par%NNZ_loc = nnz_loc
   allocate(mumps_par%IRN_loc(nnz_loc))
   allocate(mumps_par%JCN_loc(nnz_loc))
   allocate(mumps_par%A_loc(nnz_loc))
!
!..memory allocation for static condensation
   call stc_alloc
!
   if (IPRINT_TIME .eq. 1) then
      if (RANK .eq. ROOT) write(*,1003)
 1003 format(/,' STEP 2 started : Global Assembly')
      call MPI_BARRIER(mumps_par%COMM, ierr)
      start_time = MPI_Wtime()
   endif
!
!..assemble global stiffness matrix
!..loop through elements
   idec = 2
!
!$OMP PARALLEL                                  &
!$OMP PRIVATE(nrdofs,nrdofm,nrdofc,nodm,nrnodm, &
!$OMP         ndofmH,ndofmE,ndofmV,ndofmQ,      &
!$OMP         i,j,k,k1,k2,l,nod,ndof,mdle,subd)
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
!$OMP SCHEDULE(DYNAMIC)  &
!$OMP REDUCTION(+:RHS)
   do iel=1,NRELES
      mdle = mdle_list(iel)
      call get_subd(mdle, subd)
      if (subd .ne. RANK) cycle
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle,idec, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      else
         call celem(mdle,idec, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      endif
!
!  ...determine local to global dof connectivities
      l=0
!  ...H1 dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmH(i)
            l=l+1
            LCON(l) = NFIRSTH(nod)+j
         enddo
      enddo
!  ...H(curl) dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmE(i)
            l=l+1
            LCON(l) = nrdof_H + NFIRSTE(nod)+j
         enddo
      enddo
!  ...H(div) dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmV(i)
            l=l+1
            LCON(l) = nrdof_H + nrdof_E + NFIRSTV(nod)+j
         enddo
      enddo
!  ...L2 dof
      if (.not. ISTC_FLAG) then
         nod = nodm(nrnodm)
         do j=1,ndofmQ(nrnodm)
            l=l+1
            LCON(l) = nrdof_H + nrdof_E + nrdof_V + NFIRSTQ(nod)+j
         enddo
      endif
!
!  ...number of element (interface) dof
      ndof = l
!
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      do k1=1,ndof
!     ...global dof is:
         i = LCON(k1)
!     ...Assemble global load vector
         RHS(i) = RHS(i) + ZLOAD(k1)
!     ...loop through dof `to the right'
         do k2=1,ndof
!        ...global dof is:
            j = LCON(k2)
!        ...assemble
!        ...note: repeated indices are summed automatically by MUMPS
            elem_nnz_loc(iel) = elem_nnz_loc(iel) + 1_8
            k = (k1-1)*ndof + k2
            mumps_par%A_loc(  elem_nnz_loc(iel)) = ZTEMP(k)
            mumps_par%IRN_loc(elem_nnz_loc(iel)) = i
            mumps_par%JCN_loc(elem_nnz_loc(iel)) = j
         enddo
      enddo
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
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
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
      if (RANK .eq. ROOT) write(*,1009)
 1009 format(' STEP 3 started : Solve')
      call MPI_BARRIER(mumps_par%COMM, ierr)
      start_time = MPI_Wtime()
   endif
!
!..allocate RHS vector on host
   if (RANK .eq. ROOT) then
      allocate(mumps_par%RHS(mumps_par%N))
      mumps_par%RHS=ZERO
   else
      allocate(mumps_par%RHS(1))
   endif
!
!..gather RHS vector information on host
   count = mumps_par%N
   call MPI_REDUCE(RHS,mumps_par%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_par%COMM,ierr)
!
!..deallocate local RHS vector
   deallocate(RHS)
!
!..MUMPS analysis
   mumps_par%JOB = 1
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      time_stamp = MPI_Wtime()
   endif
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      time_stamp = MPI_Wtime()-time_stamp
      if (RANK .eq. ROOT) write(*,3001) time_stamp
 3001 format(' - Analysis : ',f12.5,'  seconds')
   endif
   if (mumps_par%INFOG(1) .ne. 0) then
      call mumps_destroy
      if (RANK.eq.ROOT) write(*,*) 'analysis: mumps_par%INFOG(1) .ne. 0'
      stop
   endif
   if (RANK .eq. ROOT) then
      write(*,1100) '   - MAX estimated size in GB = ',mumps_par%INFOG(16)/1000.d0
      write(*,1100) '   - SUM estimated size in GB = ',mumps_par%INFOG(17)/1000.d0
      write(*,1200) '   - Seq/parallel analysis    = ',mumps_par%INFOG(32)
      write(*,1200) '   - Ordering method used     = ',mumps_par%INFOG(7)
 1100 format(A,F11.3)
 1200 format(A,I1)
   endif
!
!..MUMPS factorization
   mumps_par%JOB = 2
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      time_stamp = MPI_Wtime()
   endif
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      time_stamp = MPI_Wtime()-time_stamp
      if (RANK .eq. ROOT) write(*,3002) time_stamp
 3002 format(' - Factorize: ',f12.5,'  seconds')
   endif
   if (mumps_par%INFOG(1) .ne. 0) then
      call mumps_destroy
      if (RANK.eq.ROOT) write(*,*) 'factorization: mumps_par%INFOG(1) .ne. 0'
      stop
   endif
   if (RANK .eq. ROOT) then
      write(*,1100) '   - MAX memory used in GB    = ',mumps_par%INFOG(21)/1000.d0
      write(*,1100) '   - SUM memory used in GB    = ',mumps_par%INFOG(22)/1000.d0
   endif
!
!..MUMPS solve
   mumps_par%JOB = 3
! 
  if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      time_stamp = MPI_Wtime()
   endif
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      time_stamp = MPI_Wtime()-time_stamp
      if (RANK .eq. ROOT) write(*,3003) time_stamp
 3003 format(' - Solve    : ',f12.5,'  seconds')
   endif
   if (mumps_par%INFOG(1) .ne. 0) then
      call mumps_destroy
      if (RANK.eq.ROOT) write(*,*) 'solve: mumps_par%INFOG(1) .ne. 0'
      stop
   endif
!
! ----------------------------------------------------------------------
!  END OF STEP 3
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
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
      if (RANK .eq. ROOT) write(*,1011)
 1011 format(' STEP 4 started : Store the solution')
      call MPI_BARRIER(mumps_par%COMM, ierr)
      start_time = MPI_Wtime()
   endif
!
!..allocate global solution vector on every processor
   allocate(SOL(mumps_par%N))
!..transfer solution from MUMPS structure to solution vector on host
   if (RANK .eq. ROOT) then
      SOL(1:mumps_par%N) = mumps_par%RHS(1:mumps_par%N)
   endif
   deallocate(mumps_par%RHS)
!..broadcast global solution from host to other processes
   count = mumps_par%N; src = ROOT
   call MPI_BCAST (SOL,count,MPI_VTYPE,src,mumps_par%COMM,ierr)
!
!$OMP PARALLEL DO                    &
!$OMP PRIVATE(i,k1,ndof,mdle,subd)   &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      mdle = mdle_list(iel)
      call get_subd(mdle, subd)
      if (subd .ne. RANK) cycle
!
      ndof = CLOC(iel)%ni
!
      allocate(ZSOL_LOC(ndof)); ZSOL_LOC=ZERO
!
      do k1=1,ndof
         i = CLOC(iel)%con(k1)
         ZSOL_LOC(k1) = SOL(i)
      enddo
      deallocate(CLOC(iel)%con)
!
      call solout(iel,ndof,NR_RHS,ZERO,ZSOL_LOC)
      deallocate(ZSOL_LOC)
!
   enddo
!$OMP END PARALLEL DO
!
! ----------------------------------------------------------------------
!  END OF STEP 4
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(4) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1012) Mtime(4)
 1012 format(' STEP 4 finished: ',f12.5,'  seconds',/)
   endif
!
   deallocate(SOL)
   deallocate(MAXDOFS)
   deallocate(NFIRSTH,NFIRSTE,NFIRSTV)
   if (.not. ISTC_FLAG) deallocate(NFIRSTQ)
   call stc_dealloc
!
!..Destroy the instance (deallocate internal data structures)
   call mumps_destroy
!
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .ge. 1)) then
      write(*,1013) sum(Mtime(1:4))
 1013 format(' par_mumps_sc FINISHED: ',f12.5,'  seconds',/)
   endif
!
!
end subroutine par_mumps_sc
