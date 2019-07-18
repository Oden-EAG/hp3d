! -----------------------------------------------------------------------
!
!    routine name       - pardiso_sc
!
! -----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - interface for PARDISO solver
!                       - routine computes global stiffness matrix
!                         and global load vector and solves the global
!                         linear system with pardiso
!                       - the assembly is computed in parallel using OMP
!                       - this routine supports both computation with or
!                         without static condensation (uses module stc)
!
!    in                 - mtype: 'H':  Hermitian (for complex)/
!                                      Symmetric (for real)
!                                      case for Lapack routines
!                                'G': General case for Lapack routines
!
! ----------------------------------------------------------------------
#include "implicit_none.h"
subroutine pardiso_sc(mtype)
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
   use stc,       only: HERM_STC, CLOC,                          &
                        stc_alloc, stc_dealloc, stc_get_nrdof
   use par_mesh,  only: DISTRIBUTED,HOST_MESH
   use mpi_param, only: RANK,ROOT
!
   implicit none
!
   character, intent(in) :: mtype
!
   character :: arg
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
   integer   :: iel,mdle,i,j,k,l,k1,k2,nod
   integer   :: inz,nz,nnz,nrdof_mdl
!
!..dummy variables
   integer :: nvoid
   VTYPE   :: zvoid
! 
!..work space for celem
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
   integer, dimension(NRELES)  :: mdle_list,m_elem_inz
   integer, allocatable        :: SIA(:), SJA(:)
!
   VTYPE, allocatable :: SA(:), RHS(:)
!
   integer*8 :: t1,t2,clock_rate,clock_max
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
   if (RANK .ne. ROOT) return
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      write(*,*) 'pardiso_sc: mesh is distributed (and not on host). returning...'
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
   if (IPRINT_TIME .eq. 1) then
      write(*,1000)
1000  format(' pardiso_sc: STARTED')
      write(*,*)
   endif
!
!..TODO multiple right-hand sides
   NR_RHS = 1
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
   mdle = 0
!..non zero elements counters
   inz = 0 ; m_elem_inz(1:NRELES) = 0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(iel) = mdle
!  ...get information from celem
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle,1, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      else
         call celem(mdle,1, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      endif
!
      nz = nrdofc
      if (HERM_STC) then
         k = (nz)*(nz+1)/2
      else
         k = nz**2
      endif
!
!  ...counting for OMP
      m_elem_inz(iel) = inz
      inz = inz + k
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
   allocate(RHS(nrdof)); RHS=ZERO
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
            j = lcon(k2)
            if (HERM_STC .and. i > j) cycle
!        ...assemble
            m_elem_inz(iel) = m_elem_inz(iel) + 1
            k = (k1-1)*ndof + k2
            SA( m_elem_inz(iel)) = ZTEMP(k)
            SIA(m_elem_inz(iel)) = i
            SJA(m_elem_inz(iel)) = j
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
   call coo2csr(SIA,SJA,SA,inz, nnz)
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
   if (HERM_STC) then
      arg = 'H'
   else
      arg = 'G'
   endif
   call pardiso_solve(SIA(1:nrdof+1),SJA(1:nnz),SA(1:nnz),arg,nnz,nrdof,NR_RHS, RHS)
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
         ZSOL_LOC(k1) = RHS(i)
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
   if (.not. ISTC_FLAG) deallocate(NFIRSTQ)
   call stc_dealloc
!
   if (IPRINT_TIME .ge. 1) then
      write(*,*)
      write(*,1013) sum(Mtime(1:5))
1013  format(' pardiso_sc FINISHED: ',f12.5,'  seconds',/)
   endif
!
!
end subroutine pardiso_sc
