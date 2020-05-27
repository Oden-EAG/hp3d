! -----------------------------------------------------------------------
!
!    routine name       - coarse_solve_pardiso
!
! -----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - routine solves the coarse system of the multigrid 
!                         solver and saves the solution vectors in an 
!                         auxiliary data structure 
!
! ----------------------------------------------------------------------
#include "implicit_none.h"
   subroutine coarse_solve_pardiso
!
   use data_structure3D, ONLY: NRNODS, NODES
   use mg_data_structure
   use assembly,         only: NR_RHS, MAXDOFM, MAXDOFS,       &
                               MAXbrickH, MAXmdlbH, NRHVAR,    &
                               MAXbrickE, MAXmdlbE, NREVAR,    &
                               MAXbrickV, MAXmdlbV, NRVVAR,    &
                               MAXbrickQ, NRQVAR,              &
                               NEXTRACT, IDBC, ZDOFD, ZERO,    &
                               ALOC, BLOC, AAUX, ZAMOD, ZBMOD, &
                               NR_PHYSA, MAXNODM
   use assembly_sc
   use control,          ONLY: ISTC_FLAG
   use stc,              ONLY: HERM_STC, CLOC,                 &
                               stc_alloc, stc_dealloc, stc_get_nrdof
   use macro_grid_info,  ONLY: ZSOL_C, NRDOF_COARSE
   use patch_info,       ONLY: CGRID_VERTICES, compute_patch_mdle
   use pardiso_data
! 
   implicit none
!
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!
!..number of local element dof for each physics variable
   integer, dimension(NR_PHYSA) :: nrdofi,nrdofb
!
!..integer counters
   integer :: nrdof_H,nrdof_E,nrdof_V
   integer :: nrdofm,nrdofc,nrnodm,nrdof,nrdof_mdl,ndof
   integer :: iel,mdle,i,j,nod,l,k1,k2,k
!
!..work space for celem
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
   integer :: m_elem_inz(NRELES)
   integer :: inz, nrpatch
!
   VTYPE, allocatable :: Rhs(:)
   VTYPE              :: zvoid
!
   real*8 :: start, finish, omp_get_wtime
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
   if (.not. ISTC_FLAG) then
      write(*,*) 'coarse_solve_pardiso: ', &
      'multigrid coarse solve requires static condensation. stop.'
      stop
   endif
!
!..use Hermitian case for static condensation
   HERM_STC = .true.
!
   NR_RHS = 1
!
   MAXDOFM = (MAXbrickH-MAXmdlbH)*NRHVAR   &
           + (MAXbrickE-MAXmdlbE)*NREVAR   &
           + (MAXbrickV-MAXmdlbV)*NRVVAR
!
! ----------------------------------------------------------------------
!  STEP 1 : 1ST LOOP THROUGH ELEMENTS, 1ST CALL TO CELEM TO GET INFO
! ----------------------------------------------------------------------
! 
!..allocate required variables for celem
   allocate(NEXTRACT(MAXDOFM))
   allocate(IDBC(MAXDOFM))
   allocate(ZDOFD(MAXDOFM,NR_RHS))
   allocate(MAXDOFS(NR_PHYSA)); MAXDOFS = 0
   MAXDOFM = 0
!
!..allocate and initialise offsets
   allocate (NFIRSTH(NRNODS)); NFIRSTH = -1
   allocate (NFIRSTE(NRNODS)); NFIRSTE = -1
   allocate (NFIRSTV(NRNODS)); NFIRSTV = -1
!
!..Step 1: determine the first dof offsets for active nodes
   nrdof_H = 0; nrdof_E = 0; nrdof_V = 0
!
!--------------------------------------------------------------------
!..first loop to celem to count active vertices
!..initialise visitation flag
   call mg_reset_visit

   nrpatch = 0
!   
!..loop through coarse grid elements   
   do iel=1,GRID(1)%nreles
!
!  ...pick up mdle node number      
      mdle = GRID(1)%mdlel(iel)
!      
!  ...get information from celem
      call celem_mg(iel,-1,mdle,1, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
!      
!  ...loop through active nodes
      do i=1,nrnodm
!
!     ...pick up the node number
         nod = nodm(i)
!         
!     ...avoid repetition
         if (NODES_MG(nod)%visit.gt.0) cycle
!
!     ...check if the node is a vertex         
         if (NODES(nod)%type .eq. 'vert') then
            nrpatch = nrpatch + 1
         endif
!
!     ...raise visitation flag
         NODES_MG(nod)%visit = 1
      enddo
!
!..end of loop through coarse grid elements   
   enddo
!
!..reset visitation flag
   call mg_reset_visit
!
   allocate(CGRID_VERTICES(nrpatch))
   GRID(1)%nrpatch = nrpatch
   nrpatch = 0
!
!--------------------------------------------------------------------
!
!..non zero elements for PARDISO
   inz = 0; m_elem_inz(1:NRELES) = 0; nrdof_mdl = 0
   do iel=1,GRID(1)%nreles
      mdle = GRID(1)%mdlel(iel)
!  ...get information from celem
      call celem_mg(iel,-1,mdle,1, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
!
      k = (nrdofc*(nrdofc+1))/2
!
!  ...counting for OMP
      m_elem_inz(iel) = inz
      inz = inz + k
!
!  ...update the maximum number of local dof
      do i=1,NR_PHYSA
         MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
      enddo
!
!  ...update the maximum number of modified element dof in the expanded mode
      MAXDOFM = max0(MAXDOFM,nrdofm)
!
!  ...loop through active nodes
      do i = 1,nrnodm
!
!     ...pick up the node number
         nod = nodm(i)
!         
!     ...avoid repetition
         if (NODES_MG(nod)%visit.ne.0) cycle
!
!     ...check if the node is a vertex         
         if (NODES(nod)%type .eq. 'vert') then
            nrpatch = nrpatch + 1 
            CGRID_VERTICES(nrpatch) = nod
            NODES_MG(nod)%visit = nrpatch
         else   
!        ...raise visitation flag
            NODES_MG(nod)%visit = -1
         endif   
      enddo
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
      call stc_get_nrdof(mdle, nrdofi,nrdofb)
      nrdof_mdl = nrdof_mdl + sum(nrdofb)
!
!..end of loop through elements
   enddo
!
   deallocate(NEXTRACT,IDBC,ZDOFD)
! 
!...total number of dof is nrdof
   nrdof = nrdof_H +  nrdof_E + nrdof_V
   NRDOF_CON = nrdof
   NRDOF_TOT = nrdof + nrdof_mdl
   NRDOF_COARSE = nrdof
!
! ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM FOR MUMPS
! ----------------------------------------------------------------------
!
   PRDS_NRHS = 1
   PRDS_N = nrdof
   PRDS_TYPE = 'H' ! hermitian matrix
   call start_pardiso
!
!..memory allocation for pardiso
   allocate(PRDS_IA(inz))
   allocate(PRDS_JA(inz))
   allocate(PRDS_A(inz))
   allocate(PRDS_RHS(nrdof))  ; PRDS_RHS = ZERO
   allocate(PRDS_XSOL(nrdof)) ; PRDS_XSOL = ZERO
   allocate(PRDS_PERM(nrdof))
!
   allocate(CLOC(GRID(1)%nreles))
!
!..assemble global stiffness matrix
!..loop through elements
!$OMP PARALLEL                                  &
!$OMP PRIVATE(nrdofs,nrdofm,nrdofc,nodm,nrnodm, &
!$OMP         ndofmH,ndofmE,ndofmV,ndofmQ,      &
!$OMP         i,j,k,k1,k2,l,nod,ndof,mdle)
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
!$OMP REDUCTION(+:PRDS_RHS)
   do iel=1,GRID(1)%nreles
      mdle = GRID(1)%mdlel(iel)
      call celem_mg(iel,-1,mdle,2, nrdofs,nrdofm,nrdofc,nodm,  &
         ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
!
!  ...determine local to global dof connectivities
      l=0
!  ...H1 dof
      do i = 1,nrnodm
         nod = nodm(i)
         do j=1,ndofmH(i)
            l=l+1
            LCON(l) = NFIRSTH(nod)+j
         enddo
      enddo
!  ...H(curl) dof
      do i = 1,nrnodm
         nod = nodm(i)
         do j=1,ndofmE(i)
            l=l+1
            LCON(l) = nrdof_H + NFIRSTE(nod)+j
         enddo
      enddo
!  ...H(div) dof
      do i = 1,nrnodm
         nod = nodm(i)
         do j=1,ndofmV(i)
            l=l+1
            LCON(l) = nrdof_H + nrdof_E + NFIRSTV(nod)+j
         enddo
      enddo
!
!  ...local number of dof is ndof
      ndof = l
!
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      do k1=1,ndof
!     ...global dof is:
         i = LCON(k1)
!     ...Assemble global load vector
         PRDS_RHS(i) = PRDS_RHS(i) + ZLOAD(k1)
!     ...loop through dof `to the right'
         do k2=1,ndof
!        ...global dof is:
            j = LCON(k2)
            if (i > j) cycle
!        ...assemble
            m_elem_inz(iel) = m_elem_inz(iel) + 1
            k = (k1-1)*ndof + k2
            PRDS_A(M_elem_inz(iel)) = ZTEMP(k)
            PRDS_IA(M_elem_inz(iel)) = i
            PRDS_JA(M_elem_inz(iel)) = j
         enddo
      enddo
!
      CLOC(iel)%ni = ndof
      allocate(CLOC(iel)%con(ndof))
      CLOC(iel)%con = LCON(1:ndof)
!  ...end of loop through elements
   enddo
!$OMP END DO
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
!----------------------------------------------------------------------
!  STEP 3: convert to compressed column format
!----------------------------------------------------------------------
!
!..convert to compressed column format
   call coo2csr(PRDS_IA,PRDS_JA,PRDS_A,inz, PRDS_NZ)
!   
!----------------------------------------------------------------------
!  STEP 4: call pardiso to solve the linear system
!----------------------------------------------------------------------
!
   PRDS_PHASE = 11 ! analysis, factorization, solve
   
   call pardiso(PRDS_PT,PRDS_MAXFCT,PRDS_MNUM,PRDS_MTYPE,PRDS_PHASE,   &
                PRDS_N,PRDS_A(1:PRDS_NZ),PRDS_IA(1:nrdof+1),           &
                PRDS_JA(1:PRDS_NZ),PRDS_PERM,PRDS_NRHS,PRDS_IPARM,     &
                PRDS_MSGLVL,PRDS_RHS,PRDS_XSOL,PRDS_ERROR)

   PRDS_PHASE = 22 ! only factorization

   call pardiso(PRDS_PT,PRDS_MAXFCT,PRDS_MNUM,PRDS_MTYPE,PRDS_PHASE,   &
                PRDS_N,PRDS_A(1:PRDS_NZ),PRDS_IA(1:nrdof+1),           &
                PRDS_JA(1:PRDS_NZ),PRDS_PERM,PRDS_NRHS,PRDS_IPARM,     &
                PRDS_MSGLVL,PRDS_RHS,PRDS_XSOL,PRDS_ERROR)
   
   PRDS_PHASE = 33 ! only solve
   PRDS_IPARM(8) = 0 ! max numbers of iterative refinement steps

   call pardiso(PRDS_PT,PRDS_MAXFCT,PRDS_MNUM,PRDS_MTYPE,PRDS_PHASE,   &
                PRDS_N,PRDS_A(1:PRDS_NZ),PRDS_IA(1:nrdof+1),           &
                PRDS_JA(1:PRDS_NZ),PRDS_PERM,PRDS_NRHS,PRDS_IPARM,     &
                PRDS_MSGLVL,PRDS_RHS,PRDS_XSOL,PRDS_ERROR)
!
! ----------------------------------------------------------------------
!  END OF STEP 4
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
!  STEP 5 : STORE SOLUTION IN THE DATA STRUCTURE
! ----------------------------------------------------------------------
!
!$OMP PARALLEL
!..allocate arrays required by solout for celem
   allocate(NEXTRACT(MAXDOFM))
   allocate(IDBC(MAXDOFM))
   allocate(ZDOFD(MAXDOFM,NR_RHS))
!
!$OMP DO                   &
!$OMP PRIVATE(i,k1,ndof,mdle)   &
!$OMP SCHEDULE(DYNAMIC)
!..loop through elements
   do iel=1,GRID(1)%nreles
!
      ndof = CLOC(iel)%ni
!
      allocate(ZSOL_LOC(ndof)) ; ZSOL_LOC=ZERO
!
      do k1=1,ndof
!     ...global dof is:
         i = CLOC(iel)%con(k1)
         ZSOL_LOC(k1) = PRDS_XSOL(i)
      enddo
!
!  ...write dof to data structure
      mdle=GRID(1)%mdlel(iel)
      call solout_mg(iel,-1,mdle,ndof,NR_RHS,ZSOL_LOC)
!
      deallocate(ZSOL_LOC)
!
!..end of loop through elements
   enddo
!$OMP END DO
   deallocate(NEXTRACT,IDBC,ZDOFD)
!$OMP END PARALLEL
   deallocate(MAXDOFS)
   deallocate(NFIRSTH,NFIRSTE,NFIRSTV)
!
!----------------------------------------------------------------------
!  END OF STEP 5
!----------------------------------------------------------------------
!
!..if multigrid is not used release memory
   if (COARSE_SOLVER .eq. NO_CSOLVE) then 
      call finalize_pardiso
   endif 
!
   call compute_patch_mdle(1)
! 
!
   end subroutine coarse_solve_pardiso
