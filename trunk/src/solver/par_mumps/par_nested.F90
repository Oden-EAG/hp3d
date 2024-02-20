!
#include "typedefs.h"
!
#if HP3D_USE_INTEL_MKL
!
! -----------------------------------------------------------------------
!
!    routine name       - par_nested
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
!                       - nested dissection solver:
!                         OpenMP MUMPS for local subdomain problem
!                         MPI MUMPS for global interface problem
!                       - this routine supports both computation with or
!                         without static condensation (uses module stc)
!                       - additionally, the local subdomain problem is
!                         condensed for each processor in the distributed
!                         mesh, followed by a global parallel solve of
!                         the condensed subdomain interface problem
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
subroutine par_nested(mtype)
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
   use control,    only: ISTC_FLAG
   use stc,        only: HERM_STC,CLOC,stc_alloc,stc_dealloc,  &
                         stc_get_nrdof,stc_bwd
   use parameters, only: ZONE,NRRHS
   use par_mumps , only: mumps_par,mumps_bub,                  &
                         mumps_start_par,mumps_destroy_par,    &
                         mumps_start_subd,mumps_destroy_subd
   use par_mesh  , only: DISTRIBUTED,HOST_MESH
   use mpi_param , only: RANK,ROOT,NUM_PROCS
   use MPI

   use mkl_spblas
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
!
!..subdomain dof counters
   integer :: nrdof_subd(0:NUM_PROCS-1)
!
!..auxiliary data structures for nested dissection
   integer, allocatable :: NOD_SUM(:),NFIRST_SUBD(:),LCON_SUBD(:)
   VTYPE  , allocatable :: ZSOL_SUBD(:),Aii(:,:),Bi(:)
   integer :: nrdof_subd_bub,nrdof_subd_con,nrdof_bub,nrdofc_bub,nrdofc_subd
   integer :: ni,nb,k_bub
   integer(8) :: nnz_bub
   integer(8) :: elem_nnz_bub(NRELES_SUBD)
!
!..sparse MKL
   integer :: mkl_stat
   integer :: nnz_ib,k_ib
   integer :: elem_nnz_ib(NRELES_SUBD)
   integer, allocatable :: Aib_row(:), Aib_col(:)
   VTYPE  , allocatable :: Aib_val(:)
   type (SPARSE_MATRIX_T) :: Aib_sparse
   type (MATRIX_DESCR)    :: Aib_descr
!
!..timer
   real(8) :: start_time,end_time,time_stamp
!
!..info (verbose output if true)
   logical :: info = .false.
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
   if ((.not. DISTRIBUTED) .or. HOST_MESH) then
      if (RANK .eq. ROOT) then
         write(*,*) 'par_nested: mesh is not distributed (or on host).'
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
1000  format(' par_nested: STARTED')
      write(*,*)
   endif
!
!..TODO multiple right-hand sides
   NR_RHS = NRRHS
!
   call mumps_start_par
   call mumps_start_subd
!
! ----------------------------------------------------------------------
!  STEP 1 : 1ST LOOP THROUGH ELEMENTS, 1ST CALL TO CELEM TO GET INFO
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1001)
 1001 format(' STEP 1 started : Get assembly info')
      start_time = MPI_Wtime()
   endif
!
!..allocate node ownership array
   allocate(NOD_OWN(NRNODS))
   allocate(NOD_SUM(NRNODS))
!
!..compute local node ownership
!  (assumes node subdomains have previously been set)
!$OMP PARALLEL DO PRIVATE(nod,subd)
   do nod=1,NRNODS
      call get_subd(nod, subd)
      if (subd .eq. RANK) then
         NOD_OWN(nod) = subd
         NOD_SUM(nod) = 1
      else
         NOD_OWN(nod) = NUM_PROCS
         NOD_SUM(nod) = 0
      endif
   enddo
!$OMP END PARALLEL DO
!
!..compute global node ownership (smallest RANK in NOD_OWN(nod) is the node owner)
!..compute whether node is subdomain bubble or subdomain interface node
!  (if NOD_SUM(nod) > 1, then node is at the interface)
   count = NRNODS
   call MPI_ALLREDUCE(MPI_IN_PLACE,NOD_OWN,count,MPI_INTEGER,MPI_MIN,mumps_par%COMM,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,NOD_SUM,count,MPI_INTEGER,MPI_SUM,mumps_par%COMM,ierr)
!
   allocate(MAXDOFS(NR_PHYSA))
   MAXDOFS = 0; MAXDOFM = 0
!
   allocate(NFIRST_DOF(NRNODS)); NFIRST_DOF = -1
   allocate(NFIRST_SUBD(NRNODS)); NFIRST_SUBD = -1
!
!..Step 1: determine the first dof offsets for active nodes
   nrdof = 0; nrdof_mdl = 0; nrdofc_subd = 0; idec = 1
   nrdof_subd_bub = 0; nrdof_subd_con = 0
!
!..non-zero counters for offsets in subdomain interior (bubble-bubble) matrix
   nnz_bub = 0_8
   elem_nnz_bub(1:NRELES_SUBD) = 0_8
!
!..non-zero counters for offsets in subdomain interface-bubble matrix 'Aib'
   nnz_ib = 0
   elem_nnz_ib(1:NRELES_SUBD) = 0
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
!  ...update the maximum number of local dof
      do i=1,NR_PHYSA
         MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
      enddo
!  ...update the maximum number of modified element dof in the expanded mode
      MAXDOFM = max0(MAXDOFM,nrdofm)
!
      nrdofc_bub = 0
!
!  ...compute offsets for nodal dof
      do i=1,nrnodm
         nod = nodm(i)
!     ...calculate total subdomain interior counter first (including duplicates)
         if (NOD_SUM(nod) .le. 1) then
            nrdofc_bub = nrdofc_bub + ndofmH(i) + ndofmE(i) + ndofmV(i)
         endif
!     ...avoid repetition within my subdomain
         if (NFIRST_SUBD(nod).ge.0) cycle
!
!     ...update the dof counter for subdomain problem
         if (NOD_SUM(nod) .le. 1) then
            NFIRST_SUBD(nod) = nrdof_subd_bub ! offset in subdomain bubble problem
            nrdof_subd_bub = nrdof_subd_bub + ndofmH(i) + ndofmE(i) + ndofmV(i)
         else
            NFIRST_SUBD(nod) = nrdof_subd_con ! offset in subdomain interface problem
            nrdof_subd_con = nrdof_subd_con + ndofmH(i) + ndofmE(i) + ndofmV(i)
         endif
!
!     ...avoid repetition within overlaps with other subdomains
         if (NOD_OWN(nod) .ne. RANK) cycle
!
!     ...store the first dof offset, and
!     ...update the dof counter
         if (NOD_SUM(nod) .gt. 1) then
            NFIRST_DOF(nod) = nrdof
            nrdof = nrdof + ndofmH(i) + ndofmE(i) + ndofmV(i) ! offset in global interface problem
         endif
      enddo
      if (.not. ISTC_FLAG) nrdof_subd_bub = nrdof_subd_bub + ndofmQ(nrnodm)
!
!  ...compute number of middle node bubble dof (nrdof_mdl)
      if (ISTC_FLAG) then
         call stc_get_nrdof(mdle, nrdofi,nrdofb)
         nrdof_mdl = nrdof_mdl + sum(nrdofb)
      endif
!
!     The following index calculations are counting duplicate indices on purpose
!
!  ...nrdofc      = number of subdomain condensed modified element dof after compression
!  ...nrdofc_subd = dof counter for local interface problem
!  ...k           = number of subdomain condensed non-zero entries in element stiffness matrix
      nrdofc_subd = nrdofc_subd + (nrdofc-nrdofc_bub)
      k = (nrdofc-nrdofc_bub)**2
!  ...nrdofc_bub = number of subdomain interior modified element dof after compression
!  ...k_bub      = number of subdomain interior non-zero entries in element stiffness matrix 'Abb'
!  ...k_ib       = number of subdomain interface-bubble non-zero entries in element stiffness matrix 'Aib'
      k_bub = nrdofc_bub**2
      k_ib  = (nrdofc-nrdofc_bub)*nrdofc_bub
!
!  ...subdomain interior counters for OpenMP (bubble-bubble matrix)
      elem_nnz_bub(iel) = nnz_bub
      nnz_bub = nnz_bub + int8(k_bub)
!  ...subdomain interior counters for OpenMP (interface-bubble matrix)
      elem_nnz_ib(iel) = nnz_ib
      nnz_ib = nnz_ib + k_ib
 1234 format(A,I10)
   enddo
!
!..compute subdomain offset
   nrdof_subd(0:NUM_PROCS-1) = 0
   nrdof_subd(RANK) = nrdof ! subdomain interface dofs from owned nodes
   count = NUM_PROCS
   call MPI_ALLREDUCE(MPI_IN_PLACE,nrdof_subd,count,MPI_INTEGER,MPI_MAX,mumps_par%COMM,ierr)
!
!..calculate prefix sum for global offsets for global interface dofs
   nrdof = 0
   do i = 0,RANK-1
      nrdof = nrdof + nrdof_subd(i)
   enddo
!$OMP PARALLEL DO
   do i = 1,NRNODS
      if ((NOD_OWN(i).eq.RANK) .and. (NOD_SUM(i).gt.1)) NFIRST_DOF(i) = NFIRST_DOF(i) + nrdof
   enddo
!$OMP END PARALLEL DO
!
   deallocate(NOD_OWN)
!
!..communicate offsets (to receive offsets for non-owned nodes within subdomain, i.e., at subdomain interface)
   count = NRNODS
   call MPI_ALLREDUCE(MPI_IN_PLACE,NFIRST_DOF,count,MPI_INTEGER,MPI_MAX,mumps_par%COMM,ierr)
!
!..calculate total number of (interface) dofs
   nrdof = 0
   do i = 0,NUM_PROCS-1
      nrdof = nrdof + nrdof_subd(i)
   enddo
!
!..compute total number of condensed element bubble dofs
   count = 1
   call MPI_ALLREDUCE(MPI_IN_PLACE,nrdof_mdl,count,MPI_INTEGER,MPI_SUM,mumps_par%COMM,ierr)
!
!..compute total number of condensed subdomain bubble dofs
   nrdof_bub = 0; count = 1;
   call MPI_ALLREDUCE(nrdof_subd_bub,nrdof_bub,count,MPI_INTEGER,MPI_SUM,mumps_par%COMM,ierr)
!
!..compute total number of non-zeros in global interface matrix
   nnz_loc = int8(nrdof_subd_con**2); nnz = 0_8; count = 1;
   call MPI_ALLREDUCE(nnz_loc,nnz,count,MPI_INTEGER8,MPI_SUM,mumps_par%COMM,ierr)
!
!..total number of (interface) dof is nrdof
   NRDOF_CON = nrdof
   NRDOF_TOT = nrdof + nrdof_bub + nrdof_mdl
!
   if (nrdof .eq. 0) then
      deallocate(MAXDOFS,NFIRST_DOF,NFIRST_SUBD,NOD_SUM)
      if (RANK .eq. ROOT) write(*,*) 'par_nested: nrdof = 0. returning.'
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
      write(*,2010) '       Number of dof   : nrdof_con      = ', NRDOF_CON,  &
                    '                         nrdof_bub      = ', nrdof_bub,  &
                    '                         nrdof_mdl      = ', nrdof_mdl,  &
                    '                         nrdof_tot      = ', NRDOF_TOT,  &
                    '       Total interf nnz: nnz            = ', nnz
   endif
   call MPI_BARRIER(mumps_par%COMM, ierr)
   if (info) then
   write(*,2011) '[', RANK, '] Local interf dof: nrdof_subd_con = ', nrdof_subd_con, &
                       '       Local bubble dof: nrdof_subd_bub = ', nrdof_subd_bub, &
                       '       Local interf nnz: nnz_loc        = ', nnz_loc,        &
                       '       Local bubble nnz: nnz_bub        = ', nnz_bub,        &
                       '       Local interf-bub: nnz_ib         = ', nnz_ib
   else
      mumps_bub%icntl(4) = 0
      mumps_par%icntl(4) = 0
   endif

2010 format(A,I12,/,A,I12,/,A,I12,/,A,I12,/,A,I12,/)
2011 format(A,I4,A,I12,/,A,I12,/,A,I12,/,A,I12,/,A,I12,/)
   call MPI_BARRIER(mumps_par%COMM, ierr)
!
!..memory allocation for local subdomain MUMPS solve
!  TODO write warning if subdomain has NO bubble dofs,
!       then skip the bubble solve with that RANK;
!       or, possibly modify mumps_par communicator.
   if (nrdof_subd_bub .eq. 0) then
      write(*,2011) '[', RANK, '] PAR_NESTED WARNING: nrdof_subd_bub = 0, no subdomain problem.'
   endif
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1003)
 1003 format(/,' STEP 2 started : Global Assembly')
      start_time = MPI_Wtime()
   endif
!
!..Set up parameters for local MUMPS solve
   mumps_bub%N = nrdof_subd_bub
   mumps_bub%NNZ = nnz_bub
   mumps_bub%NRHS = NR_RHS + nrdof_subd_con
   mumps_bub%LRHS = nrdof_subd_bub
   allocate(mumps_bub%IRN(nnz_bub))
   allocate(mumps_bub%JCN(nnz_bub))
   allocate(mumps_bub%A(nnz_bub))
   allocate(mumps_bub%RHS(mumps_bub%N * (NR_RHS+nrdof_subd_con)))
   mumps_bub%RHS=ZERO
!
   allocate(Aii(nrdof_subd_con,nrdof_subd_con)); Aii = ZERO
   allocate(Bi (nrdof_subd_con))               ; Bi = ZERO
!
!..Create sparse matrix Aib in COO format
   allocate(Aib_val(nnz_ib))
   allocate(Aib_row(nnz_ib))
   allocate(Aib_col(nnz_ib))
!
!..memory allocation for indices of condensed subdomain solution
   allocate(LCON_SUBD(nrdof_subd_con))
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
   allocate(LCON_SUBD_CON(MAXDOFM))
   allocate(ZLOAD(MAXDOFM))
   allocate(ZTEMP(MAXDOFM**2))
!
!$OMP DO SCHEDULE(DYNAMIC)
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
         if (NOD_SUM(nod) .le. 1) then ! bubble node
            do j=1,ndofmH(i)
               l=l+1
               LCON(l) = -(NFIRST_SUBD(nod)+j) ! negative index indicates bubble dof
            enddo
         else ! interface node
            do j=1,ndofmH(i)
               l=l+1
               LCON(l) = NFIRST_DOF(nod)+j ! positive index indicates interface dof
               LCON_SUBD_CON(l) = NFIRST_SUBD(nod)+j ! index for assembling subdomain matrix
               LCON_SUBD(LCON_SUBD_CON(l)) = LCON(l) ! index for retrieving local interface problem from global solution vector
            enddo
         endif
      enddo
!  ...H(curl) dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         if (NOD_SUM(nod) .le. 1) then ! bubble node
            do j=1,ndofmE(i)
               l=l+1
               LCON(l) = -(NFIRST_SUBD(nod)+ndofmH(i)+j) ! negative index indicates bubble dof
            enddo
         else ! interface node
            do j=1,ndofmE(i)
               l=l+1
               LCON(l) = NFIRST_DOF(nod)+ndofmH(i)+j ! positive index indicates interface dof
               LCON_SUBD_CON(l) = NFIRST_SUBD(nod)+ndofmH(i)+j ! index for assembling subdomain matrix
               LCON_SUBD(LCON_SUBD_CON(l)) = LCON(l) ! index for retrieving local subdomain interface problem from global solution vector
            enddo
         endif
      enddo
!  ...H(div) dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         if (NOD_SUM(nod) .le. 1) then ! bubble node
            do j=1,ndofmV(i)
               l=l+1
               LCON(l) = -(NFIRST_SUBD(nod)+ndofmH(i)+ndofmE(i)+j) ! negative index indicates bubble dof
            enddo
         else ! interface node
            do j=1,ndofmV(i)
               l=l+1
               LCON(l) = NFIRST_DOF(nod)+ndofmH(i)+ndofmE(i)+j ! positive index indicates interface dof
               LCON_SUBD_CON(l) = NFIRST_SUBD(nod)+ndofmH(i)+ndofmE(i)+j ! index for assembling subdomain matrix
               LCON_SUBD(LCON_SUBD_CON(l)) = LCON(l) ! index for retrieving local subdomain interface problem from global solution vector
            enddo
         endif
      enddo
!  ...L2 dof
      if (.not. ISTC_FLAG) then
         nod = nodm(nrnodm) ! always bubble node
         do j=1,ndofmQ(nrnodm)
            l=l+1
            LCON(l) = -(NFIRST_DOF(nod)+ndofmH(nrnodm)+ndofmE(nrnodm)+ndofmV(nrnodm)+j)
         enddo
      endif
!  ...number of element (interface) dof
      ndof = l
!
!  ...assemble the local subdomain stiffness matrices and load vector
!     (dense assembly for subdomain interface, sparse assembly for subdomain interior)
!  ...loop through element dof
!$OMP CRITICAL
      do k1=1,ndof
!     ...global interface (or subdomain bubble) dof is:
         i = LCON(k1)
         if (i .gt. 0) then ! interface dof
            !mumps_par%RHS(i) = mumps_par%RHS(i) + ZLOAD(k1)
!        ...Assemble local subdomain interface load vector
            i = LCON_SUBD_CON(k1)
            Bi(i) = Bi(i) + ZLOAD(k1)
         else ! bubble dof
!        ...Assemble local subdomain interior load vector
            mumps_bub%RHS(-i) = mumps_bub%RHS(-i) + ZLOAD(k1)
         endif
      enddo
      do k1=1,ndof ! row-major (ZTEMP)
!     ...global interface (or subdomain bubble) dof is:
         i = LCON(k1)
!     ...loop through dof `to the right'
         if (i .gt. 0) then ! row: interface dof
            do k2=1,ndof ! col
!           ...global dof is:
               j = LCON(k2)
               if (j .gt. 0) then ! col: interface dof
!              ...Assemble local subdomain interface-interface stiffness matrix 'Aii'
!              ...note: repeated indices are summed automatically by MUMPS
                  !elem_nnz_loc(iel) = elem_nnz_loc(iel) + 1_8
                  k = (k1-1)*ndof + k2
                  !mumps_par%A_loc(  elem_nnz_loc(iel)) = ZTEMP(k)
                  !mumps_par%IRN_loc(elem_nnz_loc(iel)) = i
                  !mumps_par%JCN_loc(elem_nnz_loc(iel)) = j
                  i = LCON_SUBD_CON(k1) ! subdomain interface row
                  j = LCON_SUBD_CON(k2) ! subdomain interface col
                  Aii(i,j) = Aii(i,j) + ZTEMP(k)
               else ! col: bubble dof
!              ...Assemble local subdomain interface-bubble stiffness matrix 'Aib''
                  k = (k1-1)*ndof + k2
                  i = LCON_SUBD_CON(k1) ! subdomain interface row
                                        ! -j is subdomain bubble col
                  !Aib(i,-j) = Aib(i,-j) + ZTEMP(k)
!              ...sparse assembly
                  elem_nnz_ib(iel) = elem_nnz_ib(iel) + 1
                  Aib_val(elem_nnz_ib(iel)) = ZTEMP(k)
                  Aib_row(elem_nnz_ib(iel)) = i
                  Aib_col(elem_nnz_ib(iel)) = -j
               endif
            enddo
         else ! row: bubble dof
            do k2=1,ndof
!           ...global dof is:
               j = LCON(k2)
               if (j .gt. 0) then ! col: interface dof
!              ...Assemble local subdomain bubble-interface stiffness matrix 'Abi'
                  k = (k1-1)*ndof + k2
                                        ! -i is subdomain bubble row
                  j = LCON_SUBD_CON(k2) ! subdomain interface col
                  l = mumps_bub%N*NR_RHS + (j-1)*(mumps_bub%N) - i
                  mumps_bub%RHS(l) = mumps_bub%RHS(l) + ZTEMP(k)
               else ! col: bubble dof
!              ...Assemble local subdomain bubble-bubble stiffness matrix 'Abb'
                  elem_nnz_bub(iel) = elem_nnz_bub(iel) + 1_8
                  k = (k1-1)*ndof + k2
                  ! -i is subdomain bubble row
                  ! -j is subdomain bubble col
                  mumps_bub%A(  elem_nnz_bub(iel)) = ZTEMP(k)
                  mumps_bub%IRN(elem_nnz_bub(iel)) = -i
                  mumps_bub%JCN(elem_nnz_bub(iel)) = -j
               endif
            enddo
         endif
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
   deallocate(ZBMOD,ZAMOD,LCON,LCON_SUBD_CON,ZLOAD,ZTEMP)
!$OMP END PARALLEL
!
   deallocate(MAXDOFS,NFIRST_DOF,NFIRST_SUBD,NOD_SUM)
!
   if (IPRINT_TIME .eq. 1) then
      !write(*,2101) RANK,MPI_Wtime()-start_time
 2101 format('[',I4,'] - Subdomain Assembly: ',f12.5,'  seconds')
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(2) =  end_time-start_time
      if (RANK .eq. ROOT) write(*,1004) Mtime(2)
 1004 format(' STEP 2 finished: ',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  STEP 3: call mumps to solve the local linear system
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1005)
 1005 format(' STEP 3 started : Local solve')
      start_time = MPI_Wtime()
   endif
!
!..MUMPS analysis
   mumps_bub%JOB = 1
!
   if (IPRINT_TIME .eq. 1) time_stamp = MPI_Wtime()
#if C_MODE
   call zmumps(mumps_bub)
#else
   call dmumps(mumps_bub)
#endif
   if (mumps_bub%INFO(1) .ne. 0) then
      call mumps_destroy_subd
      write(*,*) '[',RANK,'] analysis: mumps_bub%INFO(1) .ne. 0'
      stop
   endif
   if (IPRINT_TIME .eq. 1 .and. info) then
      time_stamp = MPI_Wtime()-time_stamp
      write(*,3101) RANK,time_stamp,                                       &
         '         - estimated size in GB = ',mumps_bub%INFO(15)/1000.d0,  &
         '         - ordering method used = ',mumps_bub%INFOG(7)
 3101 format('[',I4,'] - Analysis : ',f12.5,'  seconds',/,A,F11.3,/,A,I1)
   endif
!
   if (info) call MPI_BARRIER(mumps_par%COMM, ierr)
!
  35 continue
!
!..MUMPS factorization
   mumps_bub%JOB = 2
!
   if (IPRINT_TIME .eq. 1) time_stamp = MPI_Wtime()
#if C_MODE
   call zmumps(mumps_bub)
#else
   call dmumps(mumps_bub)
#endif
   if (mumps_bub%INFO(1) .ne. 0) then
      write(*,*) '[',RANK,'] factorization: mumps_bub%INFO(1) .ne. 0'
      if (mumps_bub%INFO(1) .eq. -9) then
         write(*,*) '[',RANK,'] Increasing workspace, trying factorization again...'
         mumps_bub%icntl(14) = mumps_bub%icntl(14) + 10 ! increase workspace by 10 percent
         goto 35
      endif
      call mumps_destroy_subd
      stop
   endif
   if (IPRINT_TIME .eq. 1 .and. info) then
      time_stamp = MPI_Wtime()-time_stamp
      write(*,3102) RANK,time_stamp,   &
         '         - memory used in GB    = ',mumps_bub%INFO(22)/1000.d0
 3102 format('[',I4,'] - Factorize: ',f12.5,'  seconds',/,A,F11.3)
   endif
!
   if (info) call MPI_BARRIER(mumps_par%COMM, ierr)
!
!..MUMPS solve
   mumps_bub%JOB = 3
!
   if (IPRINT_TIME .eq. 1) time_stamp = MPI_Wtime()
#if C_MODE
   call zmumps(mumps_bub)
#else
   call dmumps(mumps_bub)
#endif
   if (mumps_bub%INFO(1) .ne. 0) then
      call mumps_destroy_subd
      write(*,*) '[',RANK,'] solve: mumps_bub%INFO(1) .ne. 0'
      stop
   endif
   if (IPRINT_TIME .eq. 1 .and. info) then
      time_stamp = MPI_Wtime()-time_stamp
      write(*,3103) RANK,time_stamp
 3103 format('[',I4,'] - Solve    : ',f12.5,'  seconds')
   endif
!
   if (IPRINT_TIME .eq. 1) then
      !write(*,3104) RANK,MPI_Wtime()-start_time
 3104 format('[',I4,'] - Local Solve: ',f12.5,'  seconds')
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(3) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1006) Mtime(3)
 1006 format(' STEP 3 finished: ',f12.5,'  seconds',/)
   endif
!
   deallocate(mumps_bub%IRN,mumps_bub%JCN,mumps_bub%A)
!
!----------------------------------------------------------------------
!  STEP 4: construct condensed global system from subdomains
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1007)
 1007 format(' STEP 4 started : Construct condensed global system')
      start_time = MPI_Wtime()
   endif
!
   ni = nrdof_subd_con; nb = mumps_bub%N
!
   call MPI_BARRIER(mumps_par%COMM, ierr); time_stamp = MPI_Wtime()
!
   Aib_descr%type = SPARSE_MATRIX_TYPE_GENERAL
#if C_MODE
   mkl_stat = MKL_SPARSE_Z_CREATE_COO(Aib_sparse,              &
                                      SPARSE_INDEX_BASE_ONE,   &
                                      ni,nb,nnz_ib,            &
                                      Aib_row,Aib_col,Aib_val)
   mkl_stat = MKL_SPARSE_Z_MM(SPARSE_OPERATION_NON_TRANSPOSE,              &
                              -ZONE,                                       &
                              Aib_sparse,Aib_descr,                        &
                              SPARSE_LAYOUT_COLUMN_MAJOR,                  &
                              mumps_bub%RHS(1:nb*NR_RHS),                  &
                              NR_RHS, nb, ZONE, Bi, ni)
   mkl_stat = MKL_SPARSE_Z_MM(SPARSE_OPERATION_NON_TRANSPOSE,              &
                              -ZONE,                                       &
                              Aib_sparse,Aib_descr,                        &
                              SPARSE_LAYOUT_COLUMN_MAJOR,                  &
                              mumps_bub%RHS(nb*NR_RHS+1:nb*(NR_RHS+ni)),   &
                              ni, nb, ZONE, Aii, ni)
#else
   mkl_stat = MKL_SPARSE_D_CREATE_COO(Aib_sparse,              &
                                      SPARSE_INDEX_BASE_ONE,   &
                                      ni,nb,nnz_ib,            &
                                      Aib_row,Aib_col,Aib_val)
   mkl_stat = MKL_SPARSE_D_MM(SPARSE_OPERATION_NON_TRANSPOSE,              &
                              -ZONE,                                       &
                              Aib_sparse,Aib_descr,                        &
                              SPARSE_LAYOUT_COLUMN_MAJOR,                  &
                              mumps_bub%RHS(1:nb*NR_RHS),                  &
                              NR_RHS, nb, ZONE, Bi, ni)
   mkl_stat = MKL_SPARSE_D_MM(SPARSE_OPERATION_NON_TRANSPOSE,              &
                              -ZONE,                                       &
                              Aib_sparse,Aib_descr,                        &
                              SPARSE_LAYOUT_COLUMN_MAJOR,                  &
                              mumps_bub%RHS(nb*NR_RHS+1:nb*(NR_RHS+ni)),   &
                              ni, nb, ZONE, Aii, ni)
#endif
!
   call MPI_BARRIER(mumps_par%COMM, ierr)
   if (RANK.eq.ROOT) then
      write(*,7891) '  Sparse: ', MPI_Wtime()-time_stamp
7891  format(A,F12.5)
   endif
!
   mkl_stat = MKL_SPARSE_DESTROY(Aib_sparse)
   deallocate(Aib_val,Aib_row,Aib_col)
!
!..use 64bit parallel analysis if nnz > 2B
!  (sequential metis/scotch using 32bit currently)
   if (nnz > 2.14e9_8) mumps_par%icntl(28) = 2
!
!..percentage increase in estimated workspace for global interface problem
   mumps_par%icntl(14) = 30
!
!..memory allocation for global MUMPS solve
   mumps_par%N = NRDOF_CON
   mumps_par%NNZ_loc = nnz_loc
   mumps_par%NRHS = NR_RHS ! currently assumed one right-hand side only
   ! mumps_par%LRHS = NRDOF_CON
   allocate(mumps_par%IRN_loc(mumps_par%NNZ_loc))
   allocate(mumps_par%JCN_loc(mumps_par%NNZ_loc))
   allocate(mumps_par%A_loc(mumps_par%NNZ_loc))
   allocate(mumps_par%RHS(mumps_par%N * NR_RHS))
   mumps_par%RHS=ZERO
!
!..convert dense local interface matrix into sparse global interface matrix
!  and convert dense local load vector into dense global load vector
!$OMP PARALLEL DO PRIVATE(k2,i,j,l)
   do k1 = 1,nrdof_subd_con ! col
      i = LCON_SUBD(k1)
      mumps_par%RHS(i) = Bi(k1)
      do k2 = 1,nrdof_subd_con ! row
         j = LCON_SUBD(k2)
         l = (k1-1)*nrdof_subd_con + k2
         mumps_par%A_loc(l) = Aii(k2,k1)
         mumps_par%IRN_loc(l) = j
         mumps_par%JCN_loc(l) = i
      enddo
   enddo
!$OMP END PARALLEL DO
!
   deallocate(Aii,Bi)
!
!..gather RHS vector information on host
   count = mumps_par%N
   if (RANK .eq. ROOT) then
      call MPI_REDUCE(MPI_IN_PLACE,mumps_par%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_par%COMM,ierr)
   else
      call MPI_REDUCE(mumps_par%RHS,mumps_par%RHS,count,MPI_VTYPE,MPI_SUM,ROOT,mumps_par%COMM,ierr)
   endif
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(4) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1008) Mtime(4)
 1008 format(' STEP 4 finished: ',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  STEP 5: call mumps to solve the global linear system
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1009)
 1009 format(' STEP 5 started : Global solve')
      start_time = MPI_Wtime()
   endif
!
   call par_solve(mumps_par)
   !call par_fiber(mumps_par,nrdof_subd,NUM_PROCS,1)
!
  if (IPRINT_TIME .eq. 1) then
     call MPI_BARRIER(mumps_par%COMM, ierr)
     end_time = MPI_Wtime()
     Mtime(5) =  end_time-start_time
     if (RANK .eq. ROOT) write(*,1010) Mtime(5)
1010 format(' STEP 5 finished: ',f12.5,'  seconds',/)
  endif
!
   if (associated(mumps_par%A_loc))   deallocate(mumps_par%A_loc)
   if (associated(mumps_par%IRN_loc)) deallocate(mumps_par%IRN_loc)
   if (associated(mumps_par%JCN_loc)) deallocate(mumps_par%JCN_loc)
!
!----------------------------------------------------------------------
!  STEP 6: compute bubble solution from condensed solution dofs
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1011)
 1011 format(' STEP 6 started : Backsubstitution for bubble solution')
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
!..Construct local interface solution from global solution vector
!  uses connectivity from subdomain interface dof index to global interface dof index
   allocate(ZSOL_SUBD(nrdof_subd_con))
!$OMP PARALLEL DO PRIVATE(i)
   do k1 = 1,nrdof_subd_con
      i = LCON_SUBD(k1) ! get global index of k1-th subdomain interface dof
      ZSOL_SUBD(k1) = mumps_par%RHS(i)
   enddo
!$OMP END PARALLEL DO
   deallocate(LCON_SUBD)
!
!  NOTE: All vectors and matrices here are dense
!..1. Compute (dense)  vector: x_b =  -[A_bb^-1 A_bi] x_i + [A_bb^-1 * l_b]
!     and store the solution in mumps_bub%RHS(1:mumps_bub%N*NR_RHS)
   ni = nrdof_subd_con; nb = mumps_bub%N
   call stc_bwd(ni,nb,                                      &
                mumps_bub%RHS(nb*NR_RHS+1:nb*(NR_RHS+ni)),  &
                ZSOL_SUBD, mumps_bub%RHS(1:nb*NR_RHS))
   deallocate(ZSOL_SUBD)
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(6) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1012) Mtime(6)
 1012 format(' STEP 6 finished: ',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 4 : STORE SOLUTION IN THE DATA STRUCTURE
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      if (RANK .eq. ROOT) write(*,1013)
 1013 format(' STEP 7 started : Store the solution')
      start_time = MPI_Wtime()
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
         if (i .gt. 0) then ! interface dof
            ZSOL_LOC(k1) = mumps_par%RHS(i)
         else ! bubble dof
            ZSOL_LOC(k1) = mumps_bub%RHS(-i)
         endif
      enddo
      deallocate(CLOC(iel)%con)
      call solout(iel,ndof,NR_RHS,ZERO,ZSOL_LOC)
   enddo
!$OMP END DO
   deallocate(ZSOL_LOC)
!$OMP END PARALLEL
!
   if (IPRINT_TIME .eq. 1) then
      call MPI_BARRIER(mumps_par%COMM, ierr)
      end_time = MPI_Wtime()
      Mtime(7) = end_time-start_time
      if (RANK .eq. ROOT) write(*,1014) Mtime(7)
 1014 format(' STEP 7 finished: ',f12.5,'  seconds',/)
   endif
!
!..Deallocate element static condensation data structures
   call stc_dealloc
!
!..Destroy the instance (deallocate internal data structures)
   call mumps_destroy_par
   call mumps_destroy_subd
!
   call MPI_BARRIER(mumps_par%COMM, ierr)
   if ((RANK .eq. ROOT) .and. (IPRINT_TIME .ge. 1)) then
      write(*,1015) sum(Mtime(1:7))
 1015 format(' par_nested FINISHED: ',f12.5,'  seconds',/)
   endif
!
!
end subroutine par_nested

#else

subroutine par_nested(mtype)
   use mpi_param , only: RANK,ROOT
   implicit none
   character, intent(in) :: mtype
   if (RANK .eq. ROOT) then
      write(*,*) 'par_nested: HP3D_USE_INTEL_MKL = 0. Dependency is required.'
   endif
end subroutine par_nested

#endif
