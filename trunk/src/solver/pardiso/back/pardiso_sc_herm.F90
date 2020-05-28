! -----------------------------------------------------------------------
!
!    routine name       - pardiso_sc_herm
!
! -----------------------------------------------------------------------
!
!    latest revision    - Aug 2018
!
!    purpose            - routine computes global stiffness matrix
!                         and global load vector and solves the global
!                         linear system with pardiso after performing
!                         static condensation of the interior dofs
!                       - The assembly is computed in parallel using OMP
!
! ----------------------------------------------------------------------
subroutine pardiso_sc_herm
!
   use data_structure3D, ONLY: NRNODS, NRELES,NODES
   use assembly,         ONLY: NR_RHS, MAXDOFM, MAXbrickH,MAXbrickE,       &
                               MAXbrickV, MAXbrickQ, NRHVAR, NREVAR,       &
                               NRVVAR, NRQVAR, MAXDOFS, MAXDOFC, NEXTRACT, &
                               IDBC, ZDOFD, ZERO, BLOC, AAUX,ALOC, ZBMOD,  &
                               ZAMOD, NR_PHYSA, MAXNODM
   use control, only: ISTC_FLAG
   use stc
!
   implicit none
!
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!
!..integer counters
   integer   :: nrdof_H,  nrdof_E,  nrdof_V
   integer   :: nrdof_Hb, nrdof_Eb, nrdof_Vb, nrdof_Q
   integer   :: nrdofm, nrdofc, nrnodm, nrdof, ndof
   integer   :: iel, mdle,i,j,nod,l,k1,k2
   integer   :: inz,nz,nnz,knz,k
   integer   :: nQ,nHc,nEc,nVc,nHb,nEb,nVb,n1,n2,nH,nE,nV
!
!..work space for celem
   integer :: nodm(MAXNODM),  ndofmH(MAXNODM), &
              ndofmE(MAXNODM),ndofmV(MAXNODM),ndofmQ(MAXNODM)
   integer    M_elem_nz(NRELES), M_elem_inz(NRELES)
   integer    mdle_list(NRELES)
   integer, allocatable    :: SIA(:), SJA(:)
!
   VTYPE              :: zvoid
   VTYPE, allocatable :: SA(:), RHS(:)
!
!..number of variables for each physics attribute for an element
   integer(kind=8) :: t1,t2,clock_rate,clock_max
   real*8 :: norm_inf, tol,start, end, omp_get_wtime
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1000)
1000  format(' pardiso_sc_herm: STARTED')
      write(*,*)
   endif
!
   NR_RHS = 1
   MAXDOFM = MAXbrickH*NRHVAR + MAXbrickE*NREVAR  &
           + MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
!
! ----------------------------------------------------------------------
!  STEP 1 : 1ST LOOP THROUGH ELEMENTS, 1ST CALL TO CELEM TO GET INFO
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1001)
1001  format(' STEP 1 started : Get assembly info')
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
   nrdof_H = 0 ; nrdof_E = 0; nrdof_V = 0; nrdof_Q = 0
   nrdof_Hb = 0 ; nrdof_Eb = 0; nrdof_Vb = 0
!
   mdle = 0
!..non zero elements counters
   inz = 0 ; M_elem_inz = 0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(iel) = mdle
!  ...get information from celem
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle,1, nrdofs,nrdofm,nrdofc,nodm, &
                 ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
         write(*,*) 'pardiso for stc not yet implemented. stop.'
         stop
      else
         call celem(mdle,1, nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      endif
!
      nz = nrdofc-ndofmH(nrnodm)-ndofmE(nrnodm)-ndofmV(nrnodm)-ndofmQ(nrnodm)
!
      inz = inz + (nz)*(nz+1)/2
!
!  ...counting for OMP
      M_elem_nz(iel) = (nz)*(nz+1)/2
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
!
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
!..end of loop through elements
   enddo
!
   deallocate(NEXTRACT,IDBC,ZDOFD)
!
!..total number of dof is nrdof
   nrdof = nrdof_H +  nrdof_E + nrdof_V
   NRDOF_CON = nrdof
   NRDOF_TOT = NRDOF_CON + nrdof_Hb + nrdof_Eb + nrdof_Vb + nrdof_Q
!
!
! ----------------------------------------------------------------------
!  END OF STEP 1
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
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
   if (IPRINT_TIME .eq. 1) then
      write(*,1003)
1003  format(' STEP 2 started : Global Assembly')
      call system_clock( t1, clock_rate, clock_max )
   endif
!
!..memory allocation for assembly
   allocate(SIA(inz))
   allocate(SJA(inz))
   allocate(SA(inz))
   allocate(RHS(nrdof)) ; RHS=ZERO
!
   allocate(CLOC(NRELES))
!
!..assemble global stiffness matrix
!..loop through elements
!
!$OMP PARALLEL 									&
!$OMP PRIVATE(nrdofs,nrdofm,nrdofc,nodm,nrnodm, &
!$OMP         ndofmH,ndofmE,ndofmV,ndofmQ,lcon, &
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
   allocate(ZLOAD(MAXDOFM),ZTEMP(MAXDOFM**2))
!
!$OMP DO              	&
!$OMP SCHEDULE(DYNAMIC) &
!$OMP REDUCTION(+:RHS)
   do iel=1,NRELES
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle_list(iel),2, nrdofs,nrdofm,nrdofc,nodm, &
                    ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      else
         call celem(mdle_list(iel),2, nrdofs,nrdofm,nrdofc,nodm, &
                    ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      endif
!
!  ...determine local to global dof connectivities
      l=0; nHc = 0; nEc = 0; nVc = 0;
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
      CLOC(iel)%con = lcon
!
      allocate(ZTEMPc(n1**2))
      allocate(ZLOADC(n1))
!
      if ((.not ISTC_FLAG) .and. (n2 .ne. 0)) then
         call stc_herm(iel,ZTEMP(1:ndof**2),ZLOAD(1:ndof),ndof,n1,n2, &
                     nHc,nHb,nEc,nEb,nVc,nVb,nQ,ZTEMPc,ZLOADc)
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
         RHS(i) = RHS(i) + ZLOADc(k1)
!     ...loop through dof `to the right'
         do k2=1,n1
!        ...global dof is:
            j = lcon(k2)
            if (i.le.j) then
!           ...assemble
               M_elem_inz(iel) = M_elem_inz(iel) + 1
               k = (k1-1)*n1 + k2
               SA(M_elem_inz(iel))  = ZTEMPc(k)
               SIA(M_elem_inz(iel)) = i
               SJA(M_elem_inz(iel)) = j
            endif
         enddo
      enddo
!
      deallocate(ZTEMPc,ZLOADc)
!  ...end of loop through elements
   enddo
!$OMP END DO
!
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
   if (IPRINT_TIME .eq. 1) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(2) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1004) Mtime(2)
1004  format(' STEP 2 finished: ',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  STEP 3: convert to compressed column format
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1007)
1007  format(' STEP 3 started : Convert to CSR')
      call system_clock( t1, clock_rate, clock_max )
   endif
!
   call coo2csr(SIA(1:inz),SJA(1:inz),SA(1:inz),inz, nnz)
!
! ----------------------------------------------------------------------
!  END OF STEP 3
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(3) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1008) Mtime(3)
1008  format(' STEP 3 finished: ',f12.5,'  seconds',/)
   endif
!
!----------------------------------------------------------------------
!  STEP 4: call pardiso to solve the linear system
!----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1009)
1009  format(' STEP 4 started : Solve')
      call system_clock( t1, clock_rate, clock_max )
   endif
!
   call pardiso_solve(SIA(1:nrdof+1),SJA(1:nnz),SA(1:nnz),'H',nnz,nrdof,1,RHS(1:nrdof))
!
! ----------------------------------------------------------------------
!  END OF STEP 4
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(4) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1010) Mtime(4)
1010  format(' STEP 4 finished: ',f12.5,'  seconds',/)
   endif
!
   deallocate(SIA,SJA,SA)
!
! ----------------------------------------------------------------------
!  STEP 5 : STORE SOLUTION IN THE DATASTRACTURE
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1011)
1011  format(' STEP 5 started : Store the solution')
      call system_clock( t1, clock_rate, clock_max )
   endif
!
!$OMP PARALLEL
!..allocate arrays required by solout for celem
   allocate(NEXTRACT(MAXDOFM))
   allocate(IDBC(MAXDOFM))
   allocate(ZDOFD(MAXDOFM,NR_RHS))
!
!..loop through elements
!$OMP DO                                   &
!$OMP PRIVATE(k1,i,ndof,n1,n2,nHc,nEc,nVc, &
!$OMP         nHb,nEb,nVb,nQ,nH,nE,nV)     &
!$OMP SCHEDULE(DYNAMIC)
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
      allocate(ZSOL_LOCc(n1))  ; ZSOL_LOCc=ZERO
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
         write(*,*) 'condensed solution:'
         write(*,*) ZSOL_LOCc
         write(*,*) 'bubble solution:'
         write(*,*) ZSOL_LOCb
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
!..end of loop through elements
   enddo
!$OMP END DO
   deallocate(NEXTRACT,IDBC,ZDOFD)
!$OMP END PARALLEL
!
! ----------------------------------------------------------------------
!  END OF STEP 5
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      call system_clock( t2, clock_rate, clock_max )
      Mtime(5) =  real(t2 - t1,8)/real(clock_rate,8)
      write(*,1012) Mtime(5)
1012  format(' STEP 5 finished: ',f12.5,'  seconds',/)
   endif
!
   deallocate(RHS)
   deallocate(MAXDOFS)
   deallocate(NFIRSTV,NFIRSTH,NFIRSTE)
   deallocate(CLOC)
!
   if (IPRINT_TIME .ge. 1) then
      write(*,*)
      write(*,1013) sum(Mtime(1:5))
1013  format(' pardiso_sc_herm FINISHED: ',f12.5,'  seconds',/)
   endif
!
!
end subroutine pardiso_sc_herm
