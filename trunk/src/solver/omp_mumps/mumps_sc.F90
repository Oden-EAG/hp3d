!
#include "implicit_none.h"
! -----------------------------------------------------------------------
!
!    routine name       - mumps_sc
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 2019
!
!    purpose            - interface for OpenMP MUMPS solver
!                       - routine computes global stiffness matrix
!                         and global load vector and solves the global
!                         linear system with mumps
!                       - the assembly is computed in parallel using OMP
!                       - this routine supports both computation with or
!                         without static condensation (uses module stc)
!
!    in                 - mtype: 'H':  Hermitian (for complex)/
!                                      Symmetric (for real)
!                                      case for Lapack routines
!                                'G':  General case for Lapack routines
!
! ----------------------------------------------------------------------
subroutine mumps_sc(mtype)
!
   use data_structure3D, only: NRNODS, NRELES
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
   use mumps,     only: MUMPS_PAR, mumps_start, mumps_destroy
   use par_mesh,  only: DISTRIBUTED,HOST_MESH
   use mpi_param, only: RANK,ROOT
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
   integer    :: iel,mdle,i,j,k,l,k1,k2,nod
   integer(8) :: nnz
!
!..dummy variables
   integer :: nvoid
   VTYPE   :: zvoid
!
!..work space for celem
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
   integer    :: mdle_list(NRELES)
   integer(8) :: elem_nnz(NRELES)
!
   VTYPE, allocatable :: RHS(:)
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      if (RANK .eq. ROOT) then
         write(*,*) 'mumps_sc: mesh is distributed (and not on host).'
         write(*,*) 'calling par_mumps_sc...'
      endif
      call par_mumps_sc(mtype)
      return
   endif
   if (RANK .ne. ROOT) return
!
   select case(mtype)
      case('H')
         HERM_STC = .true.
      case default
         HERM_STC = .false.
   end select
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1000)
1000  format(' mumps_sc: STARTED')
      write(*,*)
   endif
!
!..TODO multiple right-hand sides
   NR_RHS = 1
   call mumps_start
!
!..case with static condensation
   if (ISTC_FLAG) then
      MAXDOFM = (MAXbrickH-MAXmdlbH)*NRHVAR   &
              + (MAXbrickE-MAXmdlbE)*NREVAR   &
              + (MAXbrickV-MAXmdlbV)*NRVVAR
!..no static condensation
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
      write(*,1001)
1001  format(' STEP 1 started : Get assembly info')
      start_time = MPI_Wtime()
   endif
!
!..allocate required variables for celem
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
   mdle = 0
!..matrix non-zero entries counter (and element offsets)
   nnz = 0_8 ; elem_nnz(1:NRELES) = 0_8
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(iel) = mdle
!  ...get information from celem
      if (ISTC_FLAG) then
         !write(*,*) 'celem_systemI, iel = ', iel
         call celem_systemI(iel,mdle,1, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      else
         !write(*,*) 'celem, iel = ', iel
         call celem(mdle,1, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      endif
!
      k = nrdofc**2
!
!  ...counting for OMP
      elem_nnz(iel) = nnz
      nnz = nnz + int8(k)
!
!  ...update the maximum number of local dof
      do i=1,NR_PHYSA
         MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
      enddo
!  ...update the maximum number of modified element dof in the expanded mode
      MAXDOFM = max0(MAXDOFM,nrdofm)
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
      if (ISTC_FLAG) then
         call stc_get_nrdof(mdle, nrdofi,nrdofb)
         nrdof_mdl = nrdof_mdl + sum(nrdofb)
      endif
!
!..end of loop through elements
   enddo
!
!..total number of (interface) dof is nrdof
   nrdof = nrdof_H +  nrdof_E + nrdof_V + nrdof_Q

   NRDOF_CON = nrdof
   NRDOF_TOT = nrdof + nrdof_mdl
!
   if (nrdof .eq. 0) then
      deallocate(MAXDOFS,NFIRSTH,NFIRSTE,NFIRSTV)
      if (.not. ISTC_FLAG) deallocate(NFIRSTQ)
      write(*,*) 'par_mumps_sc: nrdof = 0. returning.'
      return
   endif
!
! ----------------------------------------------------------------------
!  END OF STEP 1
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      end_time = MPI_Wtime()
      Mtime(1) = end_time-start_time
      write(*,1002) Mtime(1)
 1002 format(' STEP 1 finished: ',f12.5,'  seconds',/)
   endif 
!
! ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM 
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1003)
 1003 format(' STEP 2 started : Global Assembly')
      start_time = MPI_Wtime()
   endif
!
!..memory allocation for assembly
   allocate(RHS(nrdof)); RHS=ZERO
!
!..memory allocation for mumps
   MUMPS_PAR%N   = nrdof
   MUMPS_PAR%NNZ = nnz
!
   write(*,2010) ' Number of dof  : nrdof   = ', nrdof
   write(*,2010) ' Total non-zeros: nnz     = ', nnz
2010 format(A,I12)
!
   allocate(MUMPS_PAR%IRN(nnz))
   allocate(MUMPS_PAR%JCN(nnz))
   allocate(MUMPS_PAR%A(nnz))
!   
   call stc_alloc
!
!..assemble global stiffness matrix
!..loop through elements
!
!$OMP PARALLEL                                  &
!$OMP PRIVATE(nrdofs,nrdofm,nrdofc,nodm,nrnodm, &
!$OMP         ndofmH,ndofmE,ndofmV,ndofmQ,      &
!$OMP         i,j,k,k1,k2,l,nod,ndof)
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
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle_list(iel),2, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      else
         call celem(mdle_list(iel),2, nrdofs,nrdofm,nrdofc,nodm,  &
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
            elem_nnz(iel) = elem_nnz(iel) + 1
            k = (k1-1)*ndof + k2
            MUMPS_PAR%A(  elem_nnz(iel)) = ZTEMP(k)
            MUMPS_PAR%IRN(elem_nnz(iel)) = i
            MUMPS_PAR%JCN(elem_nnz(iel)) = j
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
      end_time = MPI_Wtime()
      Mtime(2) =  end_time-start_time
      write(*,1004) Mtime(2)
 1004 format(' STEP 2 finished: ',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  STEP 3: call mumps to solve the linear system
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1009)
 1009 format(' STEP 3 started : Solve')
      start_time = MPI_Wtime()
   endif
!
   allocate(MUMPS_PAR%RHS(MUMPS_PAR%N));
   MUMPS_PAR%RHS = RHS
   deallocate(RHS);
!
!#if C_MODE
!   MUMPS_PAR%JOB = 1
!   call zmumps(MUMPS_PAR)
!   MUMPS_PAR%JOB = 2
!   call zmumps(MUMPS_PAR)
!   MUMPS_PAR%JOB = 3
!   call zmumps(MUMPS_PAR)
!#else
!   MUMPS_PAR%JOB = 1
!   call dmumps(MUMPS_PAR)
!   MUMPS_PAR%JOB = 2
!   call dmumps(MUMPS_PAR)
!   MUMPS_PAR%JOB = 3
!   call dmumps(MUMPS_PAR)
!#endif
!
!..MUMPS analysis
   mumps_par%JOB = 1
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
   if (mumps_par%INFO(1) .ne. 0) then
      call mumps_destroy
      write(*,*) 'analysis: mumps_par%INFO(1) .ne. 0'
      stop
   endif
   write(*,1100) ' MUMPS: estimated size in GB = ',mumps_par%INFO(15)/1000.d0
 1100 format(A,F11.3)
!
!..MUMPS factorization
   mumps_par%JOB = 2
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
   if (mumps_par%INFO(1) .ne. 0) then
      call mumps_destroy
      write(*,*) 'factorization: mumps_par%INFO(1) .ne. 0'
      stop
   endif
!
!..MUMPS solve
   mumps_par%JOB = 3
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
   if (mumps_par%INFO(1) .ne. 0) then
      call mumps_destroy
      write(*,*) 'solve: mumps_par%INFO(1) .ne. 0'
      stop
   endif
!
! ----------------------------------------------------------------------
!  END OF STEP 3
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      end_time = MPI_Wtime()
      Mtime(3) = end_time-start_time
      write(*,1010) Mtime(3)
 1010 format(' STEP 3 finished: ',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 4 : STORE SOLUTION IN THE DATASTRACTURE
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1011)
 1011 format(' STEP 4 started : Store the solution')
      start_time = MPI_Wtime()
   endif
! 
!$OMP PARALLEL DO          &
!$OMP PRIVATE(i,k1,ndof)   &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES
!
      ndof = CLOC(iel)%ni
!
      allocate(ZSOL_LOC(ndof)); ZSOL_LOC=ZERO
!
      do k1=1,ndof
         i = CLOC(iel)%con(k1)
         ZSOL_LOC(k1) = MUMPS_PAR%RHS(i)
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
      end_time = MPI_Wtime()
      Mtime(4) = end_time-start_time
      write(*,1012) Mtime(4)
 1012 format(' STEP 4 finished: ',f12.5,'  seconds',/)
   endif
!   
   deallocate(MAXDOFS)
   deallocate(NFIRSTH,NFIRSTE,NFIRSTV)
   if (.not. ISTC_FLAG) deallocate(NFIRSTQ)
   call stc_dealloc
!
!..Destroy the instance (deallocate internal data structures)
   call mumps_destroy
!
   if (IPRINT_TIME .ge. 1) then
      write(*,*)
      write(*,1013) sum(Mtime(1:4))
 1013 format(' mumps_sc FINISHED: ',f12.5,'  seconds',/)
   endif
!
!
end subroutine mumps_sc
