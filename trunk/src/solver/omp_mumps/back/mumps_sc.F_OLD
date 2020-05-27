! -----------------------------------------------------------------------
!
!    routine name       - mumps_sc
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 17
!
!    purpose            - routine computes global stiffness matrix
!                         and global load vector and solves the global
!                         linear system with MUMPS after performing
!                         static condensation bubbles
!                       - The assembly is computed in parallel using OMP
!    in                 - mtype:  'H': Hermitian (for complex)/Symmetric
!                                      (for real) case for Lapack routines
!                                 'G': General case for Lapack routines
!
!
! ----------------------------------------------------------------------
!
subroutine mumps_sc(mtype)
!
   use data_structure3D, ONLY: NRNODS, NRELES
   use assembly,         ONLY: NR_RHS, MAXDOFM, MAXbrickH,MAXbrickE,       &
                               MAXbrickV, MAXbrickQ, NRHVAR, NREVAR,       &
                               NRVVAR, NRQVAR, MAXDOFS, MAXDOFC, NEXTRACT, &
                               IDBC, ZDOFD, ZERO, BLOC, AAUX,ALOC, ZBMOD,  &
                               ZAMOD, NR_PHYSA, MAXNODM
   use stc
   use mumps,            ONLY: MUMPS_PAR, mumps_start, mumps_destroy
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
   integer    M_elem_nz(NRELES), M_elem_inz(NRELES)
   integer    mdle_list(NRELES)
   integer(kind=8) :: t1,t2,clock_rate,clock_max, inz, knz
#if C_MODE
   complex*16              :: zvoid
   complex*16, allocatable :: RHS(:)
#else
   real*8                  :: zvoid
   real*8,     allocatable :: RHS(:)
#endif
!..type of the matrix ('H' for complex hermitian or real symmetric, 'G' for general)
   character*1, intent(in) :: mtype
   real*8                  :: norm_inf, tol, lamda
!
! ----------------------------------------------------------------------
! ---STEP 0 : SET FLAGS FOR MUMPS---------------------------------------
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1000)
 1000 format('mumps_sc: Started')
   endif
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
   if (IPRINT_TIME .eq. 1) then
      write(*,1001)
 1001 format('mumps_sc STEP 1: Get assembly info from celem')
      call system_clock( t1, clock_rate, clock_max )
   endif
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

   mdle = 0
!..non zero elements for MUMPS
   inz = 0 ; M_elem_inz = 0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(iel) = mdle
!  ...get information from celem
      call celem(mdle,1,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
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
!
! ----------------------------------------------------------------------
! END OF STEP 1
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(1) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1002) Mtime(1)
 1002 format('mumps_sc STEP 1: Finished in',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM FOR MUMPS
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1003)
 1003 format('mumps_sc STEP 2: Global Assembly')
      call system_clock( t1, clock_rate, clock_max )
   endif
!
!..memory allocation for mumps
   MUMPS_PAR%N = Nrdof
   MUMPS_PAR%NZ = inz

   allocate(MUMPS_PAR%IRN(MUMPS_PAR%NZ))
   allocate(MUMPS_PAR%JCN(MUMPS_PAR%NZ))
   allocate(MUMPS_PAR%A(MUMPS_PAR%NZ))
   allocate(RHS(MUMPS_PAR%N)); RHS = ZERO
!
   allocate(CLOC(NRELES))
!
!..assemble global stiffness matrix
!..loop through elements
!
!$OMP PARALLEL                                  &
!$OMP PRIVATE(nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
!$OMP         ndofmE,ndofmV,nrnodm,ndofmQ, &
!$OMP         i,j,k,k1,k2,l,nod,ndof,n1,n2,nHb, &
!$OMP         nEb,nVb,nQ,nHc,nEc,nVc)
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
   allocate(ZLOAD(MAXDOFM),ZTEMP(MAXDOFM**2))
!
!$OMP DO                &
!$OMP SCHEDULE(DYNAMIC) &
!$OMP REDUCTION(+:RHS)
   do iel=1,NRELES
      call celem(mdle_list(iel),2,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm, ZLOAD,ZTEMP)
!
!  ...determine local to global dof connectivities
      l=0 ; nHc = 0 ; nEc = 0 ; nVc = 0 ;
!  ...H1 dof
      do i = nrnodm-1,1,-1
         nod = nodm(i)
         do j=1,ndofmH(i)
            l=l+1 ; nHc = nHc + 1
            LCON(l) = NFIRSTH(nod)+j
         enddo
      enddo
!  ...H(curl) dof
      do i = nrnodm-1,1,-1
         nod = nodm(i)
         do j=1,ndofmE(i)
            l=l+1 ; nEc = nEc + 1
            LCON(l) = nrdof_H + NFIRSTE(nod)+j
         enddo
      enddo
!  ...H(div) dof
      do i = nrnodm-1,1,-1
         nod = nodm(i)
         do j=1,ndofmV(i)
            l=l+1 ; nVc = nVc + 1
            LCON(l) = nrdof_H + nrdof_E + NFIRSTV(nod)+j
         enddo
      enddo
!
!  ...local number of dof is ndof
      nHb = ndofmH(nrnodm)
      nEb = ndofmE(nrnodm)
      nVb = ndofmV(nrnodm)
      nQ  = ndofmQ(nrnodm)
!
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
      CLOC(iel)%con = LCON
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
         i = LCON(k1)
!     ...Assemble global load vector
         RHS(i) = RHS(i) + ZLOADc(k1)
!     ...loop through dof `to the right'
         do k2=1,n1
!        ...global dof is:
            j = LCON(k2)
!        ...assemble
            M_elem_inz(iel) = M_elem_inz(iel) + 1
            k = (k1-1)*n1 + k2
            MUMPS_PAR%A(M_elem_inz(iel)) = ZTEMPc(k)
            MUMPS_PAR%IRN(M_elem_inz(iel)) = i
            MUMPS_PAR%JCN(M_elem_inz(iel)) = j
         enddo
      enddo
      deallocate(ZTEMPc,ZLOADc)
!
!  ...end of loop through elements
   enddo
!$OMP END DO
   do i=1,NR_PHYSA
      deallocate(BLOC(i)%array)
      do j=1,1,NR_PHYSA
         deallocate(ALOC(i,j)%array)
      enddo
         deallocate(AAUX(i)%array)
   enddo
!
   deallocate(NEXTRACT,IDBC,ZDOFD,BLOC,AAUX,ALOC,ZBMOD,ZAMOD)
   deallocate(LCON,ZLOAD,ZTEMP)
!$OMP END PARALLEL
!
!----------------------------------------------------------------------
!  END OF STEP 2
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(2) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1004) Mtime(2)
 1004 format('mumps_sc STEP 2: Finished in',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  END OF STEP 3
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1006)
 1006 format('mumps_sc STEP 3: N/A')
   endif
!
! ----------------------------------------------------------------------
!  STEP 4 : CALL MUMPS TO SOLVE
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1007)
 1007 format('mumps_sc STEP 4: Solve')
      call system_clock( t1, clock_rate, clock_max )
   endif
!
   allocate(MUMPS_PAR%RHS(MUMPS_PAR%N));
   MUMPS_PAR%RHS = RHS
   deallocate(RHS);
!
#if C_MODE
   MUMPS_PAR%JOB = 1
   call zmumps(MUMPS_PAR)
   MUMPS_PAR%JOB = 2
   call zmumps(MUMPS_PAR)
#else
   MUMPS_PAR%JOB = 1
   call dmumps(MUMPS_PAR)
   MUMPS_PAR%JOB = 2
   call dmumps(MUMPS_PAR)
#endif
!
   ! call sparse_inverse_iteration(nrdof,lamda)
   ! write(*,*) '1: minimum lamda = ',lamda
   ! call sparse_inverse_iteration_indefinite(nrdof,lamda)
   ! write(*,*) '2: minimum lamda = ',lamda
!
#if C_MODE
   MUMPS_PAR%JOB = 3
   call zmumps(MUMPS_PAR)
#else
   MUMPS_PAR%JOB = 3
   call dmumps(MUMPS_PAR)
#endif
!
!----------------------------------------------------------------------
!  END OF STEP 4
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(4) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1008) Mtime(4)
 1008 format('mumps_sc STEP 4: Finished in',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 5 : STORE SOLUTION IN THE DATASTRACTURE
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1009)
 1009 format('mumps_sc STEP 5: Store the solution')
      call system_clock( t1, clock_rate, clock_max )
   endif
!
!$OMP PARALLEL
!
!..allocate arrays required by solout for celem
   allocate(NEXTRACT(MAXDOFM))
   allocate(IDBC(MAXDOFM))
   allocate(ZDOFD(MAXDOFM,NR_RHS))
!
!..loop through elements
!
!$OMP DO                                     &
!$OMP PRIVATE(k1,i,ndof,n1,n2,nHc,nEc,nVc,   &
!$OMP         nHb,nEb,nVb,nQ,nH,nE,nV)
   do iel=1,NRELES
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
      allocate(ZSOL_LOC(ndof)) ; ZSOL_LOC=ZERO
      allocate(ZSOL_LOCc(n1)) ; ZSOL_LOCc=ZERO
!
      do k1=1,n1
!     ...global dof is:
         i = CLOC(iel)%con(k1)
         ZSOL_LOCc(k1) = MUMPS_PAR%RHS(i)
      enddo
!
      deallocate(CLOC(iel)%con)
!
      if (n2 .ne. 0) then
         allocate(ZSOL_LOCb(n2))
         call stc_in(iel,n1,n2,ZSOL_LOCc,ZSOL_LOCb)
         ZSOL_LOC(1:nHb) = ZSOL_LOCb(1:nHb)
         ZSOL_LOC(nHb+1:nH) = ZSOL_LOCc(1:nHc)
         ZSOL_LOC(nH+1:nH+nEb) = ZSOL_LOCb(nHb+1:nHb+nEb)
         ZSOL_LOC(nH+nEb+1:nH+nE) = ZSOL_LOCc(nHc+1:nHc+nEc)
         ZSOL_LOC(nH+nE+1:nH+nE+nVb) = ZSOL_LOCb(nHb+nEb+1:nHb+nEb+nVb)
         ZSOL_LOC(nH+nE+nVb+1:nH+nE+nV)  = ZSOL_LOCc(nHc+nEc+1:nHc+nEc+nVc)
         ZSOL_LOC(nH+nE+nV+1:ndof) = ZSOL_LOCb(nHb+nEb+nVb+1:n2)
         deallocate(ZSOL_LOCb)
      else
         ZSOL_LOC = ZSOL_LOCc
      endif
!
      call solout(iel,ndof,NR_RHS,ZERO,ZSOL_LOC)
      deallocate(ZSOL_LOCc,ZSOL_LOC)
!
!  ...end of loop through elements
   enddo
!$OMP END DO
!
   deallocate(NEXTRACT,IDBC,ZDOFD)
!
!$OMP END PARALLEL
!
   deallocate(MAXDOFS)
   deallocate(NFIRSTV,NFIRSTH,NFIRSTE)
   deallocate(CLOC)
!
!----------------------------------------------------------------------
!  END OF STEP 5
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(5) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1010) Mtime(5)
 1010 format('mumps_sc STEP 5: Finished in',f12.5,'  seconds',/)
   endif
!
!..Destroy the instance (deallocate internal data structures)
!
   call mumps_destroy
!
   if (IPRINT_TIME .ge. 1) then
      write(*,1011) sum(Mtime(1:5))
 1011 format('mumps_sc: Finished in',f12.5,'  seconds',/)
   endif
!
end subroutine mumps_sc
