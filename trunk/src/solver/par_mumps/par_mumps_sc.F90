!
#include "typedefs.h"
! -----------------------------------------------------------------------
!
!    routine name       - par_mumps_sc
!
! -----------------------------------------------------------------------
!
!    latest revision    - Mar 2020
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
   use data_structure3D, only: NRNODS, NRELES_SUBD, ELEM_SUBD, &
                               get_subd
   use assembly,         only: NR_RHS, MAXDOFM, MAXDOFS,       &
                               MAXbrickH, MAXmdlbH,            &
                               MAXbrickE, MAXmdlbE,            &
                               MAXbrickV, MAXmdlbV,            &
                               MAXbrickQ,                      &
                               NEXTRACT, IDBC, ZDOFD, ZERO,    &
                               ALOC, BLOC, AAUX, ZAMOD, ZBMOD, &
                               NR_PHYSA, MAXNODM
   use assembly_sc
   use control,   only: ISTC_FLAG
   use stc,       only: HERM_STC,CLOC,                         &
                        stc_alloc,stc_dealloc,stc_get_nrdof
   use par_mumps, only: mumps_par,mumps_start_par,mumps_destroy_par
   use par_mesh , only: DISTRIBUTED,HOST_MESH
   use mpi_param, only: RANK,ROOT,NUM_PROCS
   use MPI      , only: MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,  &
                        MPI_INTEGER,MPI_INTEGER8,              &
                        MPI_REAL8,MPI_COMPLEX16
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
   integer    :: iel,mdle,subd,i,j,k,l,k1,k2,nod,idec
!
!..MPI variables
   integer :: count,src,ierr
!
!..dummy variables
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
   call mumps_start_par
!
!  ----------------------------------------------------------------------
!  STEP 1 : 1ST LOOP THROUGH ELEMENTS, 1ST CALL TO CELEM TO GET INFO
!  ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
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
   call MPI_ALLREDUCE(MPI_IN_PLACE,NOD_OWN,count,MPI_INTEGER,MPI_MIN,mumps_par%COMM,ierr)
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
!..compute offsets for owned nodes
   ! maybe OMP parallelize later
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
!  ...update the maximum number of local dof
      do i=1,NR_PHYSA
         MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
      enddo
!  ...update the maximum number of modified element dof in the expanded mode
      MAXDOFM = max0(MAXDOFM,nrdofm)
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
      if (.not. ISTC_FLAG) nrdof = nrdof + ndofmQ(nrnodm)
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
   call MPI_ALLREDUCE(MPI_IN_PLACE,nrdof_subd,count,MPI_INTEGER,MPI_MAX,mumps_par%COMM,ierr)
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
   call MPI_ALLREDUCE(MPI_IN_PLACE,NFIRST_DOF,count,MPI_INTEGER,MPI_MAX,mumps_par%COMM,ierr)
!
!..calculate total number of (interface) dofs
   nrdof = 0
   do i = 1,NUM_PROCS
      nrdof = nrdof + nrdof_subd(i)
   enddo
!
!..compute total number of condensed bubble dofs
   count = 1
   call MPI_ALLREDUCE(MPI_IN_PLACE,nrdof_mdl,count,MPI_INTEGER,MPI_SUM,mumps_par%COMM,ierr)
!
!..compute total number of non-zeros in global matrix
   count = 1
   call MPI_ALLREDUCE(nnz_loc,nnz,count,MPI_INTEGER8,MPI_SUM,mumps_par%COMM,ierr)
!
!..total number of (interface) dof is nrdof
   NRDOF_CON = nrdof
   NRDOF_TOT = nrdof + nrdof_mdl
!
   deallocate(NOD_OWN)
!
   if (nrdof .eq. 0) then
      deallocate(MAXDOFS,NFIRST_DOF)
      if (RANK .eq. ROOT) write(*,*) 'par_mumps_sc: nrdof = 0. returning.'
      return
   endif
!
!  ----------------------------------------------------------------------
!  END OF STEP 1
!  ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(1) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1002) Mtime(1)
 1002 format(' STEP 1 finished: ',f12.5,'  seconds',/)
   endif
!
!  ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM
!  ----------------------------------------------------------------------
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
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1003)
 1003 format(/,' STEP 2 started : Global Assembly')
      start_time = MPI_Wtime()
   endif
!
!..use 64bit parallel analysis if nnz > 2B
!  (sequential metis/scotch using 32bit currently)
   if (nnz > 2.14e9_8) mumps_par%icntl(28) = 2
!
!..percentage increase in estimated workspace for global interface problem
   mumps_par%icntl(14) = 30
!
!..memory allocation for MUMPS solver
   mumps_par%N = nrdof
   mumps_par%NNZ_loc = nnz_loc
   allocate(mumps_par%IRN_loc(nnz_loc))
   allocate(mumps_par%JCN_loc(nnz_loc))
   allocate(mumps_par%A_loc(nnz_loc))
   allocate(mumps_par%RHS(mumps_par%N))
   mumps_par%RHS=ZERO
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
!$OMP CRITICAL
      do k1=1,ndof
!     ...global dof is:
         i = LCON(k1)
!     ...Assemble global load vector
         mumps_par%RHS(i) = mumps_par%RHS(i) + ZLOAD(k1)
      enddo
!$OMP END CRITICAL
      do k1=1,ndof
!     ...global dof is:
         i = LCON(k1)
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
   deallocate(MAXDOFS,NFIRST_DOF)
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(2) =  end_time-start_time
      if (RANK .eq. ROOT) write(*,1004) Mtime(2)
 1004 format(' STEP 2 finished: ',f12.5,'  seconds',/)
   endif
!
!  ----------------------------------------------------------------------
!  STEP 3: CALL MUMPS TO SOLVE THE LINEAR SYSTEM
!  ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1009)
 1009 format(' STEP 3 started : Solve')
      start_time = MPI_Wtime()
   endif
!
!..gather RHS vector information on host
   count = mumps_par%N
   if (RANK .eq. ROOT) then
      call MPI_REDUCE(MPI_IN_PLACE,mumps_par%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_par%COMM,ierr)
   else
      call MPI_REDUCE(mumps_par%RHS,mumps_par%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_par%COMM,ierr)
   endif
!
   call par_solve(mumps_par)
!
   if (associated(mumps_par%A_loc))   deallocate(mumps_par%A_loc)
   if (associated(mumps_par%IRN_loc)) deallocate(mumps_par%IRN_loc)
   if (associated(mumps_par%JCN_loc)) deallocate(mumps_par%JCN_loc)
!
!  ----------------------------------------------------------------------
!  END OF STEP 3
!  ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(3) =  end_time-start_time
      if (RANK .eq. ROOT) write(*,1010) Mtime(3)
 1010 format(' STEP 3 finished: ',f12.5,'  seconds',/)
   endif
!
!  ----------------------------------------------------------------------
!  STEP 4 : STORE SOLUTION IN THE DATA STRUCTURE
!  ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1011)
 1011 format(' STEP 4 started : Store the solution')
      start_time = MPI_Wtime()
   endif
!
!..broadcast global solution from host to other processes
   count = mumps_par%N; src = ROOT
   call MPI_BCAST (mumps_par%RHS,count,MPI_VTYPE,src,mumps_par%COMM,ierr)
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      time_stamp = MPI_Wtime()-start_time
      if (RANK .eq. ROOT) write(*,3004) time_stamp
 3004 format(' - Broadcast: ',f12.5,'  seconds')
   endif
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
         ZSOL_LOC(k1) = mumps_par%RHS(i)
      enddo
      deallocate(CLOC(iel)%con)
      call solout(iel,ndof,NR_RHS,ZERO,ZSOL_LOC)
   enddo
!$OMP END DO
   deallocate(ZSOL_LOC)
!$OMP END PARALLEL
!
!  ----------------------------------------------------------------------
!  END OF STEP 4
!  ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(4) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1012) Mtime(4)
 1012 format(' STEP 4 finished: ',f12.5,'  seconds',/)
   endif
!
   call stc_dealloc
!
!..Destroy the instance (deallocate internal data structures)
   call mumps_destroy_par
!
   call MPI_BARRIER(mumps_par%COMM, ierr)
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .ge. 1)) then
      write(*,1013) sum(Mtime(1:4))
 1013 format(' par_mumps_sc FINISHED: ',f12.5,'  seconds',/)
   endif
!
!
end subroutine par_mumps_sc
