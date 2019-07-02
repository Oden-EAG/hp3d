! -----------------------------------------------------------------------
!
!    routine name       - mumps_interf
!
! -----------------------------------------------------------------------
!
!    latest revision    - Jan 16
!
!    purpose            - routine computes global stiffness matrix
!                         and global load vector and solves the global
!                         linear system with MUMPS
!
! ----------------------------------------------------------------------
!
   subroutine mumps_interf(Number_of_RHS)
!
   use control
   use data_structure3D
   use assembly
   use assembly_sc
   use stc
   use omp_lib
#include "syscom.blk"

   integer, intent(in)  :: Number_of_RHS
!
#if C_MODE
   complex*16, allocatable :: M_b(:)
   complex*16, allocatable :: M_A (:)
#else
   real*8,     allocatable :: M_b(:)
   real*8,     allocatable :: M_A (:)
#endif
   integer,    allocatable :: M_I (:), M_J(:)

!..local solution vectors
!..work space for celem
   dimension nodm(MAXNODM),  ndofmH(MAXNODM), &
             ndofmE(MAXNODM),ndofmV(MAXNODM), &
             ndofmQ(MAXNODM)
   integer   M_elem_nz(NRELES), M_elem_inz(NRELES)
   integer   mdle_list(NRELES)
   integer(kind=8) :: t1,t2,clock_rate,clock_max
!
!..number of variables for each physics attribute for an element
   dimension nrdofs(NR_PHYSA)
!
!---------------------------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1000)
 1000 format('mumps_interf: Started')
   endif
!
   MAXDOFM = MAXbrickH*NRHVAR + MAXbrickE*NREVAR  &
          + MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
!
   MAXDOFMC = MAXDOFM
!
   NR_RHS = Number_of_RHS
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1001)
 1001 format('mumps_interf STEP 1: Get assembly info from celem')
   endif
! 
   call system_clock ( t1, clock_rate, clock_max )
!   
!..allocate required variables for celem
   allocate(NEXTRACT(MAXDOFMC))
   allocate(IDBC(MAXDOFMC))
   allocate(ZDOFD(MAXDOFMC,NR_RHS))
   allocate(MAXDOFS(NR_PHYSA)); MAXDOFS = 0
   MAXDOFM = 0; MAXDOFC = 0
!
!..allocate and initialise offsets
!
   allocate (NFIRSTH(NRNODS)); NFIRSTH = -1
   allocate (NFIRSTE(NRNODS)); NFIRSTE = -1
   allocate (NFIRSTV(NRNODS)); NFIRSTV = -1
   allocate (NFIRSTQ(NRNODS)); NFIRSTQ = -1
!
!..Step 1: determine the first dof offsets for active nodes
   nrdof_H = 0 ; nrdof_E = 0; nrdof_V = 0; nrdof_Q = 0
   mdle = 0
!..non zero elements for MUMPS
   nz = 0 ; inz = 0 ; M_elem_inz = 0
!
!
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdle_list(iel) = mdle
!  ...get information from celem
      call celem(mdle,1,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm, zvoid,zvoid)
      nz = nz + nrdofc**2
!
!  ...counting for OMP
      M_elem_nz(iel) = nrdofc**2
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
      do i = nrnodm,1,-1
         if (ndofmH(i).gt.0) then
            nod = nodm(i)
!        ...avoid repetition
            if (NFIRSTH(nod).ge.0) cycle
!        ...store the first dof offset
            NFIRSTH(nod) = nrdof_H
!        ...update the H1 dof counter
            nrdof_H = nrdof_H + ndofmH(i)
         endif
      enddo
!
!  ...compute offsets for H(curl) dof
      do i = nrnodm,1,-1
         if (ndofmE(i).gt.0) then
            nod = nodm(i)
!        ...avoid repetition
            if (NFIRSTE(nod).ge.0) cycle
!        ...store the first dof offset
            NFIRSTE(nod) = nrdof_E
!        ...update the H(curl) dof counter
            nrdof_E = nrdof_E + ndofmE(i)
         endif
      enddo
!
!  ...compute offsets for H(div) dof
      do i=nrnodm,1,-1
         if (ndofmV(i).gt.0) then
            nod = nodm(i)
!        ...avoid repetition
            if (NFIRSTV(nod).ge.0) cycle
!        ...store the first dof offset
            NFIRSTV(nod) = nrdof_V
!        ...update the H(div) dof counter
            nrdof_V = nrdof_V + ndofmV(i)
         endif
      enddo
!
!  ...compute offsets for L2 dof
      i = nrnodm
      if (ndofmQ(i).gt.0) then
         nod = nodm(i)
!     ...avoid repetition
         if (NFIRSTQ(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTQ(nod) = nrdof_Q
!     ...update the L2 dof counter
         nrdof_Q = nrdof_Q + ndofmQ(i)
      endif
!
!....end of loop through elements
  enddo
!
! 
 deallocate(NEXTRACT,IDBC,ZDOFD)
! 
!...total number of dof is nrdof
  nrdof = nrdof_H +  nrdof_E + nrdof_V +  nrdof_Q
  NRDOF_TOT = nrdof
  NRDOF_CON = nrdof

!
! ----------------------------------------------------------------------
! END OF STEP 1 
! ----------------------------------------------------------------------
!
   call system_clock ( t2, clock_rate, clock_max )
   Mtime(1) =  real(t2 - t1,8)/real(clock_rate,8)

   if (IPRINT_TIME .eq. 1) then
      write(*,1002) Mtime(1)
 1002 format('mumps_interf STEP 1: Finished in',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 2 : ASSEMBLE AND STORE IN SPARSE FORM FOR MUMPS
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1003)
 1003 format('mumps_interf STEP 2: Global Assembly')
   endif
   call system_clock ( t1, clock_rate, clock_max )
!..ALLOCATE VECTORS FOR MUMPS
!
   ALLOCATE (M_A(NZ)); M_A = ZERO
   ALLOCATE (M_I(NZ)); M_I = 0
   ALLOCATE (M_J(NZ)); M_J = 0
   ALLOCATE (M_b(nrdof*NR_RHS));  M_b = ZERO
!
!..STEP 2 compute element matrices, assemble global stiffness matrix
!
!..assemble global stiffness matrix
!..loop through elements
!$OMP PARALLEL
   allocate(NEXTRACT(MAXDOFMC))
   allocate(IDBC(MAXDOFMC))
   allocate(ZDOFD(MAXDOFMC,NR_RHS))
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
      allocate(ZLOAD(MAXDOFM*NR_RHS),ZTEMP(MAXDOFM**2))
      ZLOAD=ZERO; ZTEMP=ZERO
!$OMP DO PRIVATE(nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
!$OMP            ndofmE,ndofmV,nrnodm,ndofmQ,k,k1,k2,l,i,nod,ndof,ii)
   do iel=1,NRELES
!  ...compute element matrices
      call celem(mdle_list(iel),2,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm, ZLOAD,ZTEMP)
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
      i = nrnodm
      nod = nodm(i)
      do j=1,ndofmQ(i)
         l=l+1
         LCON(l) = nrdof_H + nrdof_E + nrdof_V + NFIRSTQ(nod)+j
      enddo
!  ...local number of dof is ndof
       ndof = l
      if (ndof.gt.max_el_dof) then
         write(*,*) 'mumps_interf : WARNING : increase max_el_dof for lcon '
         write(*,8020) ndof, max_el_dof
8020     format('mumps_interf : ndof,max_el_dof =  ', i5,1x,i5)
      endif

!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      do k1=1,ndof
!     ...global dof is:
         i = LCON(k1)
!     ...Assemble global load vector
         do ii = 1, NR_RHS
            M_b(i+(ii-1)*nrdof) = M_b(i+(ii-1)*nrdof) + ZLOAD(k1+(ii-1)*ndof)
         enddo   

!     ...loop through dof `to the right'
         do k2=1,ndof
!        ...global dof is:
            j = LCON(k2)
!        ...assemble
            M_elem_inz(iel) = M_elem_inz(iel) + 1
            k = (k1-1)*ndof + k2
            M_A(M_elem_inz(iel)) = ZTEMP(k)
            M_I(M_elem_inz(iel)) = i
            M_J(M_elem_inz(iel)) = j
         enddo
      enddo
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
   call system_clock ( t2, clock_rate, clock_max )
   Mtime(2) =  real(t2 - t1,8)/real(clock_rate,8)

   if (IPRINT_TIME .eq. 1) then
      write(*,1004) Mtime(2)
 1004 format('mumps_interf STEP 2: Finished in',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 4 : CALL MUMPS TO SOLVE
! ----------------------------------------------------------------------
!
   if (IPRINT_TIME .eq. 1) then
      write(*,1007)
 1007 format('mumps_interf STEP 4: Solve')
   endif
   call system_clock ( t1, clock_rate, clock_max )
!
!..interface with MUMPS ...............................................
   call mumps_solve(nz,M_A,M_I,M_J,M_b,nrdof,NR_RHS)
!
   call system_clock ( t2, clock_rate, clock_max )
   Mtime(4) =  real(t2 - t1,8)/real(clock_rate,8)
   if (IPRINT_TIME .eq. 1) then
      write(*,1008) Mtime(4)
 1008 format('mumps_interf STEP 4: Finished in',f12.5,'  seconds',/)
   endif
!
! ----------------------------------------------------------------------
!  STEP 5 : STORE SOLUTION IN THE DATASTRACTURE
! ----------------------------------------------------------------------
!   
   if (IPRINT_TIME .eq. 1) then
      write(*,1009)
 1009 format('mumps_interf STEP 5: Store the solution')
   endif
    call system_clock ( t1, clock_rate, clock_max )
!
!..find local dofs
!..reconstruct global to local connectivities
!..loop through elements
!$OMP PARALLEL
   allocate(NEXTRACT(MAXDOFMC))
   allocate(IDBC(MAXDOFMC))
   allocate(ZDOFD(MAXDOFMC,NR_RHS))
   allocate(LCON(MAXDOFM))
!
!$OMP DO PRIVATE(nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
!$OMP        ndofmE,ndofmV,nrnodm,ndofmQ,k,k1,l,i,nod,ndof,zvoid)
   do iel=1,NRELES
!  ...get information from celem
      call celem(mdle_list(iel),1,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm, zvoid,zvoid)
      l=0
!  ...H1 dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmH(i)
            l=l+1
            LCON(l) = NFIRSTH(nod)+j
         enddo
      enddo
!
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
      i = nrnodm
      nod = nodm(i)
      do j=1,ndofmQ(i)
         l=l+1
         LCON(l) = nrdof_H + nrdof_E + nrdof_V + NFIRSTQ(nod)+j
      enddo
!  ...local number of dof is ndof
      ndof = l
!
      allocate(zsol_loc(ndof*NR_RHS)) ; zsol_loc=ZERO
      do ii=1,NR_RHS
         do k1=1,ndof
!        ...global dof is:
            i = LCON(k1)
            zsol_loc(k1+(ii-1)*ndof) = M_b(i+(ii-1)*nrdof)
         enddo
      enddo   
      call solout(iel,ndof,NR_RHS,1,zsol_loc)
      deallocate(zsol_loc)
!  ...end of loop through elements
   enddo
!$OMP END DO
!
   deallocate(NEXTRACT,IDBC,ZDOFD,LCON)
!$OMP END PARALLEL
   deallocate(MAXDOFS)
   deallocate(NFIRSTQ,NFIRSTV,NFIRSTH,NFIRSTE)
   deallocate(M_A,M_I,M_J)
   deallocate(M_b)
!
!
!----------------------------------------------------------------------
!  END OF STEP 5
!----------------------------------------------------------------------
!
   call system_clock ( t2, clock_rate, clock_max )
   Mtime(5) =  real(t2 - t1,8)/real(clock_rate,8)
   if (IPRINT_TIME .eq. 1) then
      write(*,1010) Mtime(5)
 1010 format('mumps_interf STEP 5: Finished in',f12.5,'  seconds',/)
   endif

   if (IPRINT_TIME .ge. 1) then
      write(*,1011) sum(Mtime(1:5))
 1011 format('mumps_interf: Finished in',f12.5,'  seconds',/)
   endif
!   
!
   end subroutine mumps_interf



! -----------------------------------------------------------------------
!
!    routine name       - mumps_solve
!
! -----------------------------------------------------------------------
!
!    latest revision    - Jan 16
!
!    purpose            - interface with MUMPS solver
!
! ----------------------------------------------------------------------
!
!
   subroutine mumps_solve(inz,M_A,M_I,M_J,M_b,Nrdof,NR_RHS)
!      
   implicit none
   include 'mpif.h'
   integer NR_RHS
#if C_MODE
   include 'zmumps_struc.h'
   type (ZMUMPS_STRUC) mumps_par
   complex*16 M_A(inz),M_b(nrdof*NR_RHS)
#else
   include 'dmumps_struc.h'
   type (DMUMPS_STRUC) mumps_par
   real*8 M_A(inz),M_b(nrdof*NR_RHS)
#endif
   integer IERR, I, Nrdof, inz
   integer M_I(inz),M_J(inz)
!-----------------------------------------------------------------------------
!
   call mpi_init(IERR)
!..Define a communicator for the package.
   mumps_par%COMM = MPI_COMM_WORLD
!..1 => host involved in factorization/solve phases, 0 otherwise
   mumps_par%PAR = 1
!..Initialize an instance of the package
   mumps_par%JOB = -1
!..for L U factorization (sym = 0, with working host)
   mumps_par%SYM = 0
!
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif


!..output for error messages
   mumps_par%icntl(1)  = 0
!
!..output for diagnostic/statistics/warning messages
   mumps_par%icntl(2)  = 0
!
!..output for global information
   mumps_par%icntl(3)  = 0
!
!..print level for error/warning/diagnostic messages
   mumps_par%icntl(4)  = 3
!
!..1 => element input format, 0 => assembled input format
   mumps_par%icntl(5)  = 0
!
!  0 => Approximate Minimum Degree (AMD)
!  2 => Approximate Minimum Fill (AMF)
!  4 => PORD
!  5 => METIS
!  6 => AMD w/ quasi-dense row detection
!  7 => automatic value
   mumps_par%icntl(7)  = 7
!..dense right hand side
   mumps_par%icntl(20) = 0
!..solution vector stored in the structure component RHS   
   mumps_par%icntl(21) = 0

!..Define problem on the host (processor 0)
   if(mumps_par%MYID .eq. 0) then
      mumps_par%N = Nrdof
      mumps_par%NZ = inz
      mumps_par%LRHS = nrdof
      mumps_par%NRHS = NR_RHS
!
      allocate(mumps_par%IRN(mumps_par%NZ))
      allocate(mumps_par%JCN(mumps_par%NZ))
      allocate(mumps_par%A(mumps_par%NZ))
      allocate(mumps_par%RHS(mumps_par%LRHS*mumps_par%NRHS))
      mumps_par%IRN = M_I
      mumps_par%JCN = M_J
      mumps_par%A   = M_A
      mumps_par%RHS = M_b
   endif
!..Call package for solution
   mumps_par%JOB = 6
!
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
!..Solution has been assembled on the host
   if (mumps_par%MYID .eq. 0 ) then
      M_b = mumps_par%RHS
   endif
!..Deallocate user data
   if ( mumps_par%MYID .eq. 0 ) then
      deallocate( mumps_par%IRN )
      deallocate( mumps_par%JCN )
      deallocate( mumps_par%A   )
      deallocate( mumps_par%RHS )
   endif
!..Destroy the instance (deallocate internal data structures)
   mumps_par%JOB = -2
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
   call mpi_finalize(IERR)
!
! 
  end subroutine mumps_solve
