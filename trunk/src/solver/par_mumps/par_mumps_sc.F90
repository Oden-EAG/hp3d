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
   use MPI      , only: MPI_COMM_WORLD,MPI_SUM,MPI_REAL8,MPI_COMPLEX16
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
   integer   :: nrdof_H,nrdof_E,nrdof_V,nrdof_Q
   integer   :: nrdofm,nrdofc,nrnodm,nrdof,ndof
   integer   :: iel,mdle,subd,i,j,k,l,k1,k2,nod,idec
   integer   :: inz,nz,nnz,nrdof_mdl
   integer   :: inz_loc
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
   integer, dimension(NRELES)  :: mdle_list,m_elem_inz,m_elem_inz_loc
!
!..workspace for right-hand side and solution vector
   VTYPE, allocatable :: RHS(:),SOL(:)
!
   integer(8) :: t1,t2,clock_rate,clock_max
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
   if ((.not. DISTRIBUTED) .or. HOST_MESH) then
      write(*,*) 'par_mumps_sc: mesh is not distributed (or on host). calling mumps_sc from host...'
      if (RANK .eq. ROOT) call mumps_sc(mtype)
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
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      write(*,1001)
1001  format(' STEP 1 started : Get assembly info')
      call system_clock( t1, clock_rate, clock_max )
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
   inz = 0; m_elem_inz(1:NRELES) = 0 ! global domain offset
   inz_loc = 0; m_elem_inz_loc(1:NRELES) = 0 ! local subdomain offset
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
!  ...nz = number of modified element dof after compression
!  ...k  = number of non-zeros entries in element stiffness matrix
      nz = nrdofc
      k = nz**2
!
!  ...global counters for OpenMP
      m_elem_inz(iel) = inz
      inz = inz + k
!
      if (subd .eq. RANK) then
!     ...subdomain counters for OpenMP
         m_elem_inz_loc(iel) = inz_loc
         inz_loc = inz_loc + k
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
      if (RANK .eq. ROOT) write(*,*) 'par_mumps_sc: nrdof = 0. returning.'
      return
   endif
!
!
! ----------------------------------------------------------------------
!  END OF STEP 1
! ----------------------------------------------------------------------
!
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(1) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1002) Mtime(1)
1002  format(' STEP 1 finished: ',f12.5,'  seconds',/)
   endif 
!
! ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM 
! ----------------------------------------------------------------------
!
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      write(*,1003)
1003  format(' STEP 2 started : Global Assembly')
      call system_clock( t1, clock_rate, clock_max )
   endif
!
!..memory allocation for load assembly
   allocate(RHS(nrdof)); RHS=ZERO
!
!..memory allocation for MUMPS solver
   mumps_par%N = nrdof
   mumps_par%NZ_loc = inz_loc
!
   write(*,2010) '  Number of dof  : nrdof   = ', nrdof
   write(*,2010) '  Total non-zeros: inz     = ', inz
   write(*,2010) '  Local non-zeros: inz_loc = ', inz_loc
2010 format(A,I10)
!
   allocate(mumps_par%IRN_loc(inz_loc))
   allocate(mumps_par%JCN_loc(inz_loc))
   allocate(mumps_par%A_loc(inz_loc))
!
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
            m_elem_inz_loc(iel) = m_elem_inz_loc(iel) + 1
            k = (k1-1)*ndof + k2
            mumps_par%A_loc(  m_elem_inz_loc(iel)) = ZTEMP(k)
            mumps_par%IRN_loc(m_elem_inz_loc(iel)) = i
            mumps_par%JCN_loc(m_elem_inz_loc(iel)) = j
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
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(2) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1004) Mtime(2)
1004  format(' STEP 2 finished: ',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  STEP 3: call mumps to solve the linear system
!----------------------------------------------------------------------
!
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      write(*,1009)
1009  format(' STEP 3 started : Solve')
      call system_clock( t1, clock_rate, clock_max )
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
#if C_MODE
   mumps_par%JOB = 1
   call zmumps(mumps_par)
   mumps_par%JOB = 2
   call zmumps(mumps_par)
   mumps_par%JOB = 3
   call zmumps(mumps_par)
#else
   mumps_par%JOB = 1
   call dmumps(mumps_par)
   mumps_par%JOB = 2
   call dmumps(mumps_par)
   mumps_par%JOB = 3
   call dmumps(mumps_par)
#endif
!
! ----------------------------------------------------------------------
!  END OF STEP 3
! ----------------------------------------------------------------------
!
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(3) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1010) Mtime(3)
1010  format(' STEP 3 finished: ',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 4 : STORE SOLUTION IN THE DATASTRACTURE
! ----------------------------------------------------------------------
!
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      write(*,1011)
1011  format(' STEP 4 started : Store the solution')
      call system_clock( t1, clock_rate, clock_max )
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
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .eq. 1)) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(4) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1012) Mtime(4)
1012  format(' STEP 4 finished: ',f12.5,'  seconds',/)
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
      write(*,*)
      write(*,1013) sum(Mtime(1:4))
1013  format(' par_mumps_sc FINISHED: ',f12.5,'  seconds',/)
   endif
!
!
end subroutine par_mumps_sc
