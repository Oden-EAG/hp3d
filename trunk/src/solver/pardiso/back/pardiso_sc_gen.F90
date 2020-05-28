! -----------------------------------------------------------------------
!
!    routine name       - pardiso_sc_gen
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 17
!
!    purpose            - routine computes global stiffness matrix
!                         and global load vector and solves the global
!                         linear system with pardiso after performing
!                         static condensation of the L2 dofs
!                       - The assembly is computed in parallel using OMP
!
! ----------------------------------------------------------------------

   subroutine pardiso_sc_gen
!
!
   use data_structure3D, ONLY: NRNODS, NRELES
   use assembly,         ONLY: NR_RHS, MAXDOFM, MAXbrickH,MAXbrickE,       &
                               MAXbrickV, MAXbrickQ, NRHVAR, NREVAR,       &
                               NRVVAR, NRQVAR, MAXDOFS, MAXDOFC, NEXTRACT, &
                               IDBC, ZDOFD, ZERO, BLOC, AAUX,ALOC, ZBMOD,  &
                               ZAMOD, NR_PHYSA, MAXNODM
   use stc
!
   implicit none
!
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!
!..integer counters
   integer :: nrdof_H,  nrdof_E,  nrdof_V
   integer :: nrdof_Hb, nrdof_Eb, nrdof_Vb, nrdof_Q
   integer :: nrdofm, nrdofc, nrnodm,nrdof, ndof
   integer :: iel, mdle,i,j,nod,l,k1,k2,k, inz, nz, nnz,knz
   integer :: nQ,nHc,nEc,nVc,nHb,nEb,nVb,n1,n2,nH,nE,nV
   real*8  :: norm_inf, tol
!
!..work space for celem
   integer :: nodm(MAXNODM),  ndofmH(MAXNODM), &
              ndofmE(MAXNODM),ndofmV(MAXNODM),ndofmQ(MAXNODM)
   integer    M_elem_nz(NRELES), M_elem_inz(NRELES)
   integer    mdle_list(NRELES)
   integer, allocatable    :: SIA(:), SJA(:)
#if C_MODE
   complex*16              :: zvoid
   complex*16, allocatable :: SA(:),  RHS(:),diag(:)
#else
   real*8                  :: zvoid
   real*8, allocatable     :: SA(:),  RHS(:),diag(:)
#endif
!
!..number of variables for each physics attribute for an element
   integer(kind=8) :: t1,t2,clock_rate,clock_max
!
!..printing flag
   if (IPRINT_TIME .eq. 1) then
      write(*,1000)
 1000 format(' pardiso_sc_gen: Started')
   endif
!
   NR_RHS = 1
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
 1001 format(' pardiso_sc_gen STEP 1: Get assembly info from celem')
   endif
!
   call system_clock ( t1, clock_rate, clock_max )
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
   nrdof_H = 0 ; nrdof_E = 0; nrdof_V = 0; nrdof_Q = 0
   nrdof_Hb = 0 ; nrdof_Eb = 0; nrdof_Vb = 0

   mdle = 0
!..non zero elements for PARDISO
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
 !    ...store the first dof offset
         NFIRSTH(nod) = nrdof_H
!     ...update the H1 dof counter
         nrdof_H = nrdof_H + ndofmH(i)
      enddo
!  ...compute offsets for H1 bubble dof
      i = nrnodm
      nod = nodm(i)
!  ...store the first dof offset
      NFIRSTH(nod) = nrdof_Hb
!  ...update the H1 dof counter
      nrdof_Hb = nrdof_Hb + ndofmH(i)
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

!  ...compute offsets for H(curl) bubble dof
      i = nrnodm
      nod = nodm(i)
!  ...store the first dof offset
      NFIRSTE(nod) = nrdof_Eb
!  ...update the H1 dof counter
      nrdof_Eb = nrdof_Eb + ndofmE(i)
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


!  ...compute offsets for L2 dof
      i = nrnodm
      nod = nodm(i)
!  ...repetition not possible
!  ...update the L2 dof counter
      nrdof_Q = nrdof_Q + ndofmQ(i)
!
!  ...end of loop through elements
   enddo
!
   deallocate(NEXTRACT,IDBC,ZDOFD)
!
!...total number of dof is nrdof
   nrdof = nrdof_H +  nrdof_E + nrdof_V
   nrdof_con = nrdof
   nrdof_tot = nrdof_con + nrdof_Hb + nrdof_Eb + nrdof_Vb + nrdof_Q
! ----------------------------------------------------------------------
!  END OF STEP 1
! ----------------------------------------------------------------------
!
   call system_clock ( t2, clock_rate, clock_max )
   Mtime(1) =  real(t2 - t1,8)/real(clock_rate,8)

   if (IPRINT_TIME .eq. 1) then
      write(*,1002) Mtime(1)
 1002 format(' pardiso_sc_gen STEP 1: Finished in',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM FOR MUMPS
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1003)
 1003 format(' pardiso_sc_gen STEP 2: Global Assembly')
   endif
   call system_clock ( t1, clock_rate, clock_max )
!
!..memory allocation for assembly
   allocate(SIA(inz))
   allocate(SJA(inz))
   allocate(SA(inz))
   allocate(RHS(nrdof)) ; RHS=ZERO
   allocate(diag(nrdof)); diag = ZERO
!
   allocate(CLOC(NRELES))
!
!..assemble global stiffness matrix
!..loop through elements
!$OMP PARALLEL
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
!$OMP DO PRIVATE(nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
!$OMP            ndofmE,ndofmV,nrnodm,ndofmQ, lcon , &
!$OMP            k,k1,k2,l,i,j,nod,ndof,n1,n2,nHb,nEb,nVb,nQ,nHc,nEc,nVc)
   do iel=1,NRELES
      call celem(mdle_list(iel),2,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm, ZLOAD,ZTEMP)
!
!  ...determine local to global dof connectivities
      l=0; nHc = 0 ; nEc = 0 ; nVc = 0 ;
!  ...H1 dof
      do i = nrnodm-1,1,-1
         nod = nodm(i)
         do j=1,ndofmH(i)
            l=l+1; nHc = nHc + 1
            lcon(l) = NFIRSTH(nod)+j
         enddo
      enddo
!  ...H(curl) dof
      do i = nrnodm-1,1,-1
         nod = nodm(i)
         do j=1,ndofmE(i)
            l=l+1; nEc = nEc + 1
            lcon(l) = nrdof_H + NFIRSTE(nod)+j
         enddo
      enddo
!  ...H(div) dof
      do i = nrnodm-1,1,-1
         nod = nodm(i)
         do j=1,ndofmV(i)
            l=l+1; nVc = nVc + 1
            lcon(l) = nrdof_H + nrdof_E + NFIRSTV(nod)+j
         enddo
      enddo
!  ...local number of dof is ndof

      nHb = ndofmH(nrnodm)
      nEb = ndofmE(nrnodm)
      nVb = ndofmV(nrnodm)
      nQ  = ndofmQ(nrnodm)

      n1 = l
      n2 = nHb + nEb + nVb + nQ
      ndof = n1 + n2

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
            call stc_gen(iel,ZTEMP(1:ndof**2),ZLOAD(1:ndof),ndof,n1,n2,  &
                 nHc,nHb,nEc,nEb,nVc,nVb,nQ,ZTEMPc,ZLOADc)
      else
         ZTEMPc(1:n1**2) = ZTEMP(1:n1**2)
         ZLOADc(1:n1) = ZLOAD(1:n1)
      endif
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      do k1=1,n1
!     ...global dof is:
         i = lcon(k1)
!     ...Assemble global load vector
!$OMP CRITICAL
         RHS(i) = RHS(i) + ZLOADc(k1)
!$OMP END CRITICAL
!     ...loop through dof `to the right'
         do k2=1,n1
!        ...global dof is:
            j = lcon(k2)
!        ...assemble
            M_elem_inz(iel) = M_elem_inz(iel) + 1
            k = (k1-1)*n1 + k2
            SA(M_elem_inz(iel))  = ZTEMPc(k)
            SIA(M_elem_inz(iel)) = i
            SJA(M_elem_inz(iel)) = j
            if (i.eq.j) then
!$OMP CRITICAL
               diag(i) = diag(i)+ZTEMPc(k)
!$OMP END CRITICAL
            endif
         enddo
      enddo
      deallocate(ZTEMPc,ZLOADc)

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
   deallocate(ZLOAD,ZTEMP)
!$OMP END PARALLEL
!

   call system_clock ( t2, clock_rate, clock_max )
   Mtime(2) =  real(t2 - t1,8)/real(clock_rate,8)

   if (IPRINT_TIME .eq. 1) then
      write(*,1004) Mtime(2)
 1004 format(' pardiso_sc_gen STEP 2: Finished in',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  STEP 3 : DIAGONAL SCALING
!----------------------------------------------------------------------
!
 !   if (IPRINT_TIME .eq. 1) then
 !      write(*,1005)
 ! 1005 format(' pardiso_sc_gen STEP 3: Diagonal scaling and zero filtering')
 !   endif

   call system_clock ( t1, clock_rate, clock_max )
!    do k = 1, inz
! !  ...get the row number
!       i = SIA(k)
! !  ...get the column number
!       j = SJA(k)
! !  ...scale by the ith diagonal element
!       SA(k)   = SA(k)/sqrt(diag(i)*diag(j))
!    enddo
!    do i = 1, nrdof
!       RHS(i) = RHS(i)/sqrt(diag(i))
!    enddo

!..filter out additional zeros. To be changed in the future
!..get the infinity norm of the matrix
   ! norm_inf = maxval(abs(SA))
   ! tol = norm_inf*1e-14
   ! knz = 0
   ! do k = 1, inz
   !    if (abs(SA(k)).ge. tol) then
   !       knz = knz+1
   !       SA(knz) = SA(k)
   !       SIA(knz) = SIA(k)
   !       SJA(knz) = SJA(k)
   !    endif
   ! enddo
   ! inz = knz
!
! ----------------------------------------------------------------------
!  END OF STEP 3
! ----------------------------------------------------------------------
!
   ! call system_clock ( t2, clock_rate, clock_max )
   Mtime(3) =  real(t2 - t1,8)/real(clock_rate,8)
 !   if (IPRINT_TIME .eq. 1) then
 !      write(*,1006) Mtime(3)
 ! 1006 format(' pardiso_sc_gen STEP 3: Finished in',f12.5,'  seconds',/)
 !   endif
!
!----------------------------------------------------------------------
!  STEP 4: convert to compressed column format
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1007)
 1007 format(' pardiso_sc_gen STEP 4: Convert to CSR')
   endif

   call system_clock ( t1, clock_rate, clock_max )

   call coo2csr(SIA(1:inz),SJA(1:inz),SA(1:inz),inz, nnz)
!
! ----------------------------------------------------------------------
!  END OF STEP 4
! ----------------------------------------------------------------------
!
   call system_clock ( t2, clock_rate, clock_max )
   Mtime(4) =  real(t2 - t1,8)/real(clock_rate,8)
   if (IPRINT_TIME .eq. 1) then
      write(*,1008) Mtime(4)
 1008 format(' pardiso_sc_gen STEP 4: Finished in',f12.5,'  seconds',/)
      write(*,2008) inz, nnz
 2008 format(' pardiso_sc_gen STEP 4: inz, nnz =', 2i16,/)

   endif
!
!----------------------------------------------------------------------
!  STEP 5: call pardiso to solve the linear system
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1009)
 1009 format(' pardiso_sc_gen STEP 5: Solve')
   endif
   call system_clock ( t1, clock_rate, clock_max )
!
   call pardiso_solve(SIA(1:nrdof+1),SJA(1:nnz),SA(1:nnz),'G',nnz,nrdof,1,RHS(1:nrdof))
!
   ! do i = 1, nrdof
   !    RHS(i) = RHS(i)/sqrt(diag(i))
   ! enddo
!
   deallocate(diag)
!
! ----------------------------------------------------------------------
!  END OF STEP 5
! ----------------------------------------------------------------------
!
   call system_clock ( t2, clock_rate, clock_max )
   Mtime(5) =  real(t2 - t1,8)/real(clock_rate,8)
   if (IPRINT_TIME .eq. 1) then
      write(*,1010) Mtime(5)
 1010 format(' pardiso_sc_gen STEP 5: Finished in',f12.5,'  seconds',/)
   endif
!
   deallocate(SIA,SJA,SA)
!
! ----------------------------------------------------------------------
!  STEP 6 : STORE SOLUTION IN THE DATASTRACTURE
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1011)
 1011 format(' pardiso_sc_gen STEP 6: Store the solution')
   endif
   call system_clock ( t1, clock_rate, clock_max )
!
!$OMP PARALLEL
!..allocate arrays required by solout for celem
   allocate(NEXTRACT(MAXDOFM))
   allocate(IDBC(MAXDOFM))
   allocate(ZDOFD(MAXDOFM,NR_RHS))
!
!$OMP DO PRIVATE(k1,i,ndof,n1,n2,nHc,nEc,nVc,nHb,nEb,nVb,nQ,nH,nE,nV)
! ...loop through elements
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
      allocate(ZSOL_LOC(ndof)); ZSOL_LOC=ZERO
      allocate(ZSOL_LOCc(n1)) ; ZSOL_LOCc=ZERO
!
      do k1=1,n1
!     ...global dof is:
         i = CLOC(iel)%con(k1)
         ZSOL_LOCc(k1) = RHS(i)
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
   deallocate(NEXTRACT,IDBC,ZDOFD)
!$OMP END PARALLEL
!
! ----------------------------------------------------------------------
!  END OF STEP 6
! ----------------------------------------------------------------------
!
   call system_clock ( t2, clock_rate, clock_max )
   Mtime(6) =  real(t2 - t1,8)/real(clock_rate,8)
   if (IPRINT_TIME .eq. 1) then
      write(*,1012) Mtime(6)
 1012 format(' pardiso_sc_gen STEP 6: Finished in',f12.5,'  seconds',/)
   endif
!
   deallocate(RHS)
   deallocate(MAXDOFS)
   deallocate(NFIRSTV,NFIRSTH,NFIRSTE)
   deallocate(CLOC)
!
   if (IPRINT_TIME .ge. 1) then
      write(*,1013) sum(Mtime(1:6))
 1013 format(' pardiso_sc_gen: Finished in',f12.5,'  seconds',/)
   endif
!
   end subroutine pardiso_sc_gen
