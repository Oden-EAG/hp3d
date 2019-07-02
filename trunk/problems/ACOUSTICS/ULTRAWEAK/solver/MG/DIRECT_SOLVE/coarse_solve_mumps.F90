! -----------------------------------------------------------------------
!
!    routine name       - coarse_solve_mumps
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 17
!
!    purpose            - routine solves the coarse system of the multigrid 
!                         solver and saves the solution vectors in an 
!                         auxiliary data structure (ZSOL_C)
!
!    in                 - mtype:  'H': Hermitian (for complex)/Symmetric
!                                      (for real) case for Lapack routines
!                                 'G': General case for Lapack routines
!
! ----------------------------------------------------------------------
!
   subroutine coarse_solve_mumps(mtype)
!      
   use data_structure3D, ONLY: NRNODS, NODES
   use mg_data_structure
   use assembly,         ONLY: NR_RHS, MAXDOFM, MAXbrickH,MAXbrickE,       &
                               MAXbrickV, MAXbrickQ, NRHVAR, NREVAR,       &
                               NRVVAR, NRQVAR, MAXDOFS, MAXDOFC, NEXTRACT, &
                               IDBC, ZDOFD, ZERO, BLOC, AAUX,ALOC, ZBMOD,  &
                               ZAMOD, NR_PHYSA, MAXNODM
   use m_assembly,       ONLY: NFIRSTH, NFIRSTE, NFIRSTV, NRDOF_CON, &
                               NRDOF_TOT, MTIME, CLOC, ZLOAD, ZLOADc, LCON,   &
                               ZTEMP, ZTEMPc, ZSOL_LOC, ZSOL_LOCc, ZSOL_LOCb
                               
   use mumps,            ONLY: MUMPS_PAR, mumps_start
   use macro_grid_info,  ONLY: ZSOL_C, NRELES_COARSE, MDLE_MACRO, NRDOF_COARSE
   use patch_info,       ONLY: NRPATCH, CGRID_VERTICES, compute_patch_mdle
! 
   implicit none
!
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!..integer counters
   integer :: nrdof_H,  nrdof_E,  nrdof_V
   integer :: nrdof_Hb, nrdof_Eb, nrdof_Vb, nrdof_Q
   integer :: nrdofm, nrdofc, nrnodm,nrdof, ndof
   integer :: iel, mdle,i,j,nod,l,k1,k2,k
   integer :: nQ,nHc,nEc,nVc,nHb,nEb,nVb,n1,n2,nH,nE,nV

!..work space for celem
   integer :: nodm(MAXNODM),  ndofmH(MAXNODM), &
              ndofmE(MAXNODM),ndofmV(MAXNODM),ndofmQ(MAXNODM)
   integer    M_elem_nz(NRELES_COARSE), M_elem_inz(NRELES_COARSE)
   integer :: inz, knz
#if C_MODE
   complex*16, allocatable :: Rhs(:)
   complex*16              :: zvoid(1)
#else
   real*8, allocatable     :: Rhs(:)
   real*8                  :: zvoid(1)
#endif
!..type of the matrix 'H' for complex hermitian or real symmetric, 'G' for general)
   character*1, intent(in) :: mtype
   real*8 :: start, finish, omp_get_wtime
!
!
! ----------------------------------------------------------------------
! ---STEP 0 : SET FLAGS FOR MUMPS---------------------------------------
! ----------------------------------------------------------------------
!
   NR_RHS = 1
   call mumps_start
!
   MAXDOFM = MAXbrickH*NRHVAR + MAXbrickE*NREVAR  &
           + MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
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
   MAXDOFM = 0; MAXDOFC = 0
!
!..allocate and initialise offsets
!
   allocate (NFIRSTH(NRNODS)); NFIRSTH = -1
   allocate (NFIRSTE(NRNODS)); NFIRSTE = -1
   allocate (NFIRSTV(NRNODS)); NFIRSTV = -1
!
!..Step 1: determine the first dof offsets for active nodes
   nrdof_H  = 0 ; nrdof_E  = 0; nrdof_V  = 0
   nrdof_Hb = 0 ; nrdof_Eb = 0; nrdof_Vb = 0; nrdof_Q = 0
!

!--------------------------------------------------------------------
!..first loop to celem to count active vertices
!..initialise visitation flag
   call mg_reset_visit

   NRPATCH = 0
!   
!..loop through coarse grid elements   
   do iel=1,NRELES_COARSE
!
!  ...pick up mdle node number      
      mdle = MDLE_MACRO(iel)
!      
!  ...get information from celem
      call celem1(mdle,1,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
!      
!  ...loop through active nodes
      do i = 1,nrnodm-1
!
!     ...pick up the node number
         nod = nodm(i)
!         
!     ...avoid repetition
         if (NODES_MG(nod)%visit.gt.0) cycle
!
!     ...check if the node is a vertex         
         if (NODES(nod)%type .eq. 'vert') then
            NRPATCH = NRPATCH + 1 
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
   allocate(CGRID_VERTICES(NRPATCH))
   NRPATCH = 0
!--------------------------------------------------------------------
! 
!..non zero elements for MUMPS
   inz = 0 ; M_elem_inz = 0
   do iel=1,NRELES_COARSE
      mdle = MDLE_MACRO(iel)
!  ...get information from celem
      call celem1(mdle,1,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
!
      nQ = ndofmQ(nrnodm)
      inz = inz + (nrdofc-ndofmH(nrnodm)-ndofmE(nrnodm)-ndofmV(nrnodm)-nQ)**2
!
!  ...counting for OMP
      M_elem_nz(iel) = (nrdofc-ndofmH(nrnodm)-ndofmE(nrnodm)-ndofmV(nrnodm)-nQ)**2
      M_elem_inz(iel) = sum(M_elem_nz(1:iel-1))
!  ...update the maximum number of local dof
      do i=1,NR_PHYSA
         MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
      enddo
!  ...update the maximum number of modified element dof in the expanded mode
      MAXDOFM = max0(MAXDOFM,nrdofm)
!  ...update the maximum number of modified element dof after compression
      MAXDOFC = max0(MAXDOFC,nrdofc)
!

!  ...loop through active nodes
      do i = 1,nrnodm-1
!
!     ...pick up the node number
         nod = nodm(i)
!         
!     ...avoid repetition
         if (NODES_MG(nod)%visit.ne.0) cycle
!
!     ...check if the node is a vertex         
         if (NODES(nod)%type .eq. 'vert') then
            NRPATCH = NRPATCH + 1 
            CGRID_VERTICES(NRPATCH) = nod
            NODES_MG(nod)%visit = NRPATCH
         else   
!        ...raise visitation flag
            NODES_MG(nod)%visit = -1
         endif   
      enddo
!
!  ...compute offsets for H1 dof
      do i = 1,nrnodm-1
         nod = nodm(i)
!     ...avoid repetition
         if (NFIRSTH(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTH(nod) = nrdof_H
!     ...update the H1 dof counter
         nrdof_H = nrdof_H + ndofmH(i)
      enddo
!
!  ...compute offsets for H1 bubble dof
      i = nrnodm
      nod = nodm(i)
!  ...store the first dof offset
      NFIRSTH(nod) = nrdof_Hb
!  ...update the H1 dof counter
      nrdof_Hb = nrdof_Hb + ndofmH(i)
!
!
!  ...compute offsets for H(curl) dof
      do i = 1,nrnodm-1
         nod = nodm(i)
!     ...avoid repetition
         if (NFIRSTE(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTE(nod) = nrdof_E
!     ...update the H(curl) dof counter
         nrdof_E = nrdof_E + ndofmE(i)
      enddo
!
!  ...compute offsets for H(curl) bubble dof
      i = nrnodm
      nod = nodm(i)
!  ...store the first dof offset
      NFIRSTE(nod) = nrdof_Eb
!  ...update the H1 dof counter
      nrdof_Eb = nrdof_Eb + ndofmE(i)
!
!
!  ...compute offsets for H(div) dof
      do i=1,nrnodm-1
         nod = nodm(i)
!     ...avoid repetition
         if (NFIRSTV(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTV(nod) = nrdof_V
!     ...update the H(div) dof counter
         nrdof_V = nrdof_V + ndofmV(i)
      enddo
!
!  ...compute offsets for H(div) bubble dof
      i = nrnodm
      nod = nodm(i)
!  ...store the first dof offset
      NFIRSTV(nod) = nrdof_Vb
!  ...update the H1 dof counter
      nrdof_Vb = nrdof_Vb + ndofmV(i)
!
!
!  ...compute offsets for L2 bubble dof
      i = nrnodm
      nod = nodm(i)
!  ...update the L2 dof counter (repetition not possible)
      nrdof_Q = nrdof_Q + ndofmQ(i)
!
!  ...end of loop through elements
   enddo
!
   deallocate(NEXTRACT,IDBC,ZDOFD)
! 
!...total number of dof is nrdof
   nrdof = nrdof_H +  nrdof_E + nrdof_V
   NRDOF_CON = nrdof
   NRDOF_TOT = NRDOF_CON + nrdof_Hb + nrdof_Eb + nrdof_Vb + nrdof_Q
   NRDOF_COARSE = NRDOF_CON

!
! ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM FOR MUMPS
! ----------------------------------------------------------------------
!
!..memory allocation for mumps
   MUMPS_PAR%N = nrdof
   MUMPS_PAR%NZ = inz

   allocate(MUMPS_PAR%IRN(MUMPS_PAR%NZ))
   allocate(MUMPS_PAR%JCN(MUMPS_PAR%NZ))
   allocate(MUMPS_PAR%A(MUMPS_PAR%NZ))
   allocate(MUMPS_PAR%RHS(MUMPS_PAR%N)) 
   allocate(Rhs(MUMPS_PAR%N)) ; Rhs = ZERO 
!
   allocate(CLOC(NRELES_COARSE))
!
!..assemble global stiffness matrix
!..loop through elements
!$omp parallel default(shared)  &
!$omp private(nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
!$omp         ndofmE,ndofmV,nrnodm,ndofmQ, lcon , &
!$omp         k,k1,k2,l,i,j,nod,ndof,n1,n2,nHb,nEb,nVb,nQ,nHc,nEc,nVc)
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
   allocate(ZLOAD(MAXDOFM),ZTEMP(MAXDOFM**2))
!$omp do reduction(+:Rhs) schedule(guided)
   do iel=1,NRELES_COARSE

      call celem1(MDLE_MACRO(iel),2,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm, ZLOAD,ZTEMP)
!
!  ...determine local to global dof connectivities
      l=0 ; nHc = 0 ; nEc = 0 ; nVc = 0 ;
!  ...H1 dof
      do i = 1,nrnodm-1
         nod = nodm(i)
         do j=1,ndofmH(i)
            l=l+1 ; nHc = nHc + 1
            lcon(l) = NFIRSTH(nod)+j
         enddo
      enddo
!  ...H(curl) dof
      do i = 1,nrnodm-1
         nod = nodm(i)
         do j=1,ndofmE(i)
            l=l+1 ; nEc = nEc + 1
            lcon(l) = nrdof_H + NFIRSTE(nod)+j
         enddo
      enddo
!  ...H(div) dof
      do i = 1,nrnodm-1
         nod = nodm(i)
         do j=1,ndofmV(i)
            l=l+1 ; nVc = nVc + 1
            lcon(l) = nrdof_H + nrdof_E + NFIRSTV(nod)+j
         enddo
      enddo
!
!  ...local number of dof is ndof

      nHb = ndofmH(nrnodm)
      nEb = ndofmE(nrnodm)
      nVb = ndofmV(nrnodm)
      nQ  = ndofmQ(nrnodm)

      n1 = l
      n2 = nHb + nEb + nVb + nQ
      ndof = n1 + n2
!
      CLOC(iel)%nHc = nHc
      CLOC(iel)%nEc = nEc
      CLOC(iel)%nVc = nVc
      CLOC(iel)%nHb = nHb
      CLOC(iel)%nEb = nEb
      CLOC(iel)%nVb = nVb
      CLOC(iel)%nQ  = nQ
!
      allocate(CLOC(iel)%con(n1))
      CLOC(iel)%con = lcon
!
      allocate(ZTEMPc(n1**2))
      allocate(ZLOADc(n1))
      if (n2 .ne. 0) then
         select case(mtype)
         case('H')
            call stc_herm(iel,ZTEMP(1:ndof**2),ZLOAD(1:ndof),ndof,n1,n2,  &
                 nHc,nHb,nEc,nEb,nVc,nVb,nQ,ZTEMPc,ZLOADc)
         case default   
            call stc_gen(iel,ZTEMP(1:ndof**2),ZLOAD(1:ndof),ndof,n1,n2,  &
                 nHc,nHb,nEc,nEb,nVc,nVb,nQ,ZTEMPc,ZLOADc)
         end select   
      else
         ZTEMPc(1:n1**2) = ZTEMP(1:n1**2)
         ZLOADc(1:n1) = ZLOAD(1:n1)
      endif
!
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      do k1=1,n1
!     ...global dof is:
         i = lcon(k1)
!     ...Assemble global load vector
         Rhs(i) = Rhs(i) + ZLOADc(k1)
!     ...loop through dof `to the right'
         do k2=1,n1
!        ...global dof is:
            j = lcon(k2)
!        ...assemble
            M_elem_inz(iel) = M_elem_inz(iel) + 1
            k = (k1-1)*n1 + k2
            MUMPS_PAR%A(M_elem_inz(iel)) = ZTEMPc(k)
            MUMPS_PAR%IRN(M_elem_inz(iel)) = i
            MUMPS_PAR%JCN(M_elem_inz(iel)) = j
         enddo
      enddo
      deallocate(ZTEMPc,ZLOADc)

!  ...end of loop through elements
   enddo
!$omp end do
   do i=1,NR_PHYSA
      deallocate(BLOC(i)%array)
      do j=1,1,NR_PHYSA
         deallocate(ALOC(i,j)%array)
      enddo
         deallocate(AAUX(i)%array)
   enddo
!
   deallocate(NEXTRACT,IDBC,ZDOFD,BLOC,AAUX,ALOC,ZBMOD,ZAMOD)
   deallocate(ZLOAD,ZTEMP)
!$omp end parallel

   MUMPS_PAR%RHS = Rhs
   deallocate(Rhs)
!
! ----------------------------------------------------------------------
!  STEP 4 : CALL MUMPS TO SOLVE
! ----------------------------------------------------------------------
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
!  STEP 5 : STORE SOLUTION IN THE DATASTRACTURE
! ----------------------------------------------------------------------
   allocate(ZSOL_C(NRELES_COARSE)) ! for prolongation
!   
!$omp parallel default(shared) &
!$omp private(k1,i,ndof,n1,n2,nHc,nEc,nVc,nHb,nEb,nVb,nQ,nH,nE,nV)
!..allocate arrays required by solout for celem
   allocate(NEXTRACT(MAXDOFM))
   allocate(IDBC(MAXDOFM))
   allocate(ZDOFD(MAXDOFM,NR_RHS))
!
!$omp do schedule(guided)
! ...loop through elements
   do iel=1,NRELES_COARSE
!
      nHc =  CLOC(iel)%nHc
      nEc =  CLOC(iel)%nEc
      nVc =  CLOC(iel)%nVc
!
      nHb =  CLOC(iel)%nHb
      nEb =  CLOC(iel)%nEb
      nVb =  CLOC(iel)%nVb
      nQ =   CLOC(iel)%nQ
!
      nH = nHc+nHb
      nE = nEc+nEb
      nV = nVc+nVb
!
      n1 = nHc + nEc + nVc
      n2 = nHb + nEb + nVb + nQ
      ndof = n1 + n2
!
      allocate(ZSOL_LOCc(n1)) ; ZSOL_LOCc=ZERO
!
      do k1=1,n1
!     ...global dof is:
         i = CLOC(iel)%con(k1)
         ZSOL_LOCc(k1) = MUMPS_PAR%RHS(i)
      enddo
!
      call store_elem_sol(iel, ZSOL_LOCc, n1)
!
      ! deallocate(CLOC(iel)%con)
      deallocate(ZSOL_LOCc)
!
!  ...end of loop through elements
   enddo
!$omp end do
   deallocate(NEXTRACT,IDBC,ZDOFD)
!$omp end parallel
   deallocate(MAXDOFS)
   deallocate(NFIRSTV,NFIRSTH,NFIRSTE)
!
!----------------------------------------------------------------------
!  END OF STEP 5
!----------------------------------------------------------------------
!
   write(*,*) 'Compute patch elements...'
   start = omp_get_wtime()

   call compute_patch_mdle
   finish = omp_get_wtime()

   write(*,1000) finish-start
1000 format(' Compute patch elements:t=         ', es13.4)    
! 
!
   end subroutine coarse_solve_mumps
