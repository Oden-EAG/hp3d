#include "typedefs.h"

module nl_solver_module
! use data_structure3D, only: NRELES
! use assembly, only: super_array
   use assembly_sc,only:NRDOF_CON
   use mumps,     only: MUMPS_PAR,mumps_destroy
! integer,save :: NUPDATE,MAXITER
! integer,save :: MODE_ASSEMBLY
real*8, save :: LINESEARCH_FACTOR,LINESEARCH_TOL !,RES_TOL 
logical,save :: LINESEARCH_FLAG, BFGS_FLAG
logical,save :: FIRST_NLSOLVE = .true.

real*8,allocatable :: RES_GLOBAL(:),RES_PREV(:),RES_CURR(:),V_BFGS(:,:),W_BFGS(:,:)!,NEW_SOL

integer(8),allocatable :: ELEM_NNZ(:),LCON0(:)
!$OMP THREADPRIVATE(LCON0)

! !..element local residual vectors
!    type(super_array), allocatable :: RLOC(:)
! !$OMP THREADPRIVATE (RLOC)

! !..modified element residual vector
!    VTYPE, allocatable :: ZRMOD(:,:)
! !$OMP THREADPRIVATE (ZRMOD)







contains

subroutine set_nl_solver_params(Lst_tmp,Lsf_tmp,Bfgs_tmp)
implicit none
real*8 , intent(in) :: Lst_tmp
logical, intent(in) :: Lsf_tmp,Bfgs_tmp


LINESEARCH_TOL = Lst_tmp
LINESEARCH_FLAG= Lsf_tmp
BFGS_FLAG      = Bfgs_tmp

end subroutine




subroutine prep_mumps(Mtype)
   use data_structure3D, only: NRNODS, NRELES, ELEM_ORDER
   use assembly,         only: NR_RHS, MAXDOFM, MAXDOFS, NR_PHYSA, MAXNODM
   use control,   only: ISTC_FLAG
   use stc,       only: HERM_STC,CLOC,stc_alloc,stc_get_nrdof
   use mumps,     only: MUMPS_PAR,mumps_start
   use assembly_sc
implicit none


!
   character, intent(in) :: mtype
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!
!..number of local element dof for each physics variable
   integer, dimension(NR_PHYSA) :: nrdofi,nrdofb
!
!..integer counters
   integer    :: nrdofm,nrdofc,nrnodm,nrdof,nrdof_mdl,ndof
   integer    :: iel,mdle,i,j,k,l,k1,k2,nod
!
!..dummy variables
   VTYPE   :: zvoid
!
!..work space for celem
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
!
!..64 bit non-zero entry counters
   integer(8) :: nnz
! 
!
   select case(Mtype)
      case('H')
         HERM_STC = .true.
      case default
         HERM_STC = .false.
   end select
!
!    if (IPRINT_TIME .eq. 1) then
!       write(*,1000)
! 1000  format(' mumps_sc: STARTED')
!       write(*,*)
!    endif
!
!..TODO multiple right-hand sides
   NR_RHS = 1
   call mumps_start
!
! ----------------------------------------------------------------------
!  STEP 1 : 1ST LOOP THROUGH ELEMENTS, 1ST CALL TO CELEM TO GET INFO
! ----------------------------------------------------------------------
!
!    if (IPRINT_TIME .eq. 1) then
!       write(*,1001)
! 1001  format(' STEP 1 started : Get assembly info')
!       start_time = MPI_Wtime()
!    endif
!
!..allocate required variables for celem
   allocate(MAXDOFS(NR_PHYSA))
   MAXDOFS = 0; MAXDOFM = 0
!
!..allocate and initialize offsets
   allocate(NFIRST_DOF(NRNODS)); NFIRST_DOF = -1
!
!..Step 1: determine the first dof offsets for active nodes
   nrdof = 0; nrdof_mdl = 0
!
!..matrix non-zero entries counter (and element offsets)
   allocate(ELEM_NNZ(NRELES))
   nnz = 0_8 ; ELEM_NNZ(1:NRELES) = 0_8
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
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
      ELEM_NNZ(iel) = nnz
      nnz = nnz + int8(k)
!
!  ...update the maximum number of local dof
      do i=1,NR_PHYSA
         MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
      enddo
!  ...update the maximum number of modified element dof in the expanded mode
      MAXDOFM = max0(MAXDOFM,nrdofm)
!
!  ...compute offsets for nodal dof
      do i=1,nrnodm
         nod = nodm(i)
!     ...avoid repetition
         if (NFIRST_DOF(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRST_DOF(nod) = nrdof
!     ...update the H(div) dof counter
         nrdof = nrdof + ndofmH(i) + ndofmE(i) + ndofmV(i)
      enddo
      if (.not. ISTC_FLAG) nrdof = nrdof + ndofmQ(nrnodm)
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
   NRDOF_CON = nrdof
   NRDOF_TOT = nrdof + nrdof_mdl
!
   if (nrdof .eq. 0) then
      deallocate(MAXDOFS,NFIRST_DOF)
      write(*,*) 'par_mumps_sc: nrdof = 0. returning.'
      return
   endif
!
! ----------------------------------------------------------------------
!  END OF STEP 1
! ----------------------------------------------------------------------
!

!
   write(*,2010) ' Number of dof  : nrdof_con = ', NRDOF_CON
   write(*,2010) '                  nrdof_tot = ', NRDOF_TOT
   write(*,2010) ' Total non-zeros: nnz       = ', nnz
2010 format(A,I12)
!
!..memory allocation for mumps
   MUMPS_PAR%N   = nrdof
   MUMPS_PAR%NNZ = nnz
   allocate(MUMPS_PAR%IRN(nnz))
   allocate(MUMPS_PAR%JCN(nnz))
   allocate(MUMPS_PAR%A(nnz))
   allocate(MUMPS_PAR%RHS(nrdof));

   allocate(LCON0(MAXDOFM))
!
!..memory allocation for static condensation
   call stc_alloc

   ! allocate(NEXTRACT(MAXDOFM))
   ! allocate(IDBC(MAXDOFM))
   ! allocate(ZDOFD(MAXDOFM,NR_RHS)) 
! $OMP PARALLEL                                  & 
! $OMP PRIVATE(nrdofs,nrdofm,nrdofc,nodm,nrnodm, &
! $OMP         ndofmH,ndofmE,ndofmV,ndofmQ,      &
! $OMP         i,j,k,k1,k2,l,mdle,nod,ndof)      &
! $OMP DO                                        &
! $OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle,1, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      else
         call celem(mdle,1, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
      endif
!
!  ...determine local to global dof connectivities
      l=0
!  ...H1 dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmH(i)
            l=l+1
            LCON0(l) = NFIRST_DOF(nod)+j
         enddo
      enddo
!  ...H(curl) dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmE(i)
            l=l+1
            LCON0(l) = NFIRST_DOF(nod)+ndofmH(i)+j
         enddo
      enddo
!  ...H(div) dof
      do i = nrnodm,1,-1
         nod = nodm(i)
         do j=1,ndofmV(i)
            l=l+1
            LCON0(l) = NFIRST_DOF(nod)+ndofmH(i)+ndofmE(i)+j
         enddo
      enddo
!  ...L2 dof
      if (.not. ISTC_FLAG) then
         nod = nodm(nrnodm)
         do j=1,ndofmQ(nrnodm)
            l=l+1
            LCON0(l) = NFIRST_DOF(nod)+ndofmH(nrnodm)+ndofmE(nrnodm)+ndofmV(nrnodm)+j
         enddo
      endif

      ! write(*,*) 'iel=',iel
      ! write(*,*) LCON0

      ndof = l
!
      CLOC(iel)%ni = ndof
      allocate(CLOC(iel)%con(ndof))
      CLOC(iel)%con = LCON0(1:ndof)
!..end of loop through elements
   enddo  
! $OMP END DO

!    deallocate(NEXTRACT,IDBC,ZDOFD)
! $OMP END PARALLEL

   deallocate(LCON0)


 !   if (IPRINT_TIME .eq. 1) then
 !      end_time = MPI_Wtime()
 !      Mtime(1) = end_time-start_time
 !      write(*,1002) Mtime(1)
 ! 1002 format(' STEP 1 finished: ',f12.5,'  seconds',/)
 !   endif


end subroutine




subroutine assembly_mumps(Isel,Neq,RHS)
! Isel =-1 --> retrieve stored OLD solution
! Isel = 0 --> retrieve stored NEW solution
! Isel = 1 --> Only assemble global residual RES_GLOBAL
! Isel = 2 --> Assemble both RES_GLOBAL and matrix A (which is stored in the mumps structure)
use data_structure3D, only: NRNODS, NRELES, ELEM_ORDER
use assembly,         only: NR_RHS, MAXDOFM, MAXDOFS,       &
                            NEXTRACT, IDBC, ZDOFD, ZERO,    &
                            ALOC, BLOC, AAUX, ZAMOD, ZBMOD, &
                            NR_PHYSA, MAXNODM
use control,   only: ISTC_FLAG
use stc,       only: HERM_STC,CLOC,stc_alloc,stc_dealloc,stc_get_nrdof
! use mumps,     only: 
use assembly_sc
implicit none
integer,intent(in) :: Isel,Neq
VTYPE  ,intent(out):: RHS(Neq)

!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!
!..number of local element dof for each physics variable
   integer, dimension(NR_PHYSA) :: nrdofi,nrdofb
!
!..integer counters
   integer    :: nrdofm,nrdofc,nrnodm,nrdof,nrdof_mdl,ndof
   integer    :: iel,mdle,i,j,k,l,k1,k2,nod
!
!..dummy variables
   VTYPE   :: zvoid
!
!..work space for celem
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
!
! 
   RHS = ZERO

   if (IPRINT_TIME .eq. 1) then
      ! write(*,1003)
 ! 1003 format(/,' STEP 2 started : Global Assembly')
      select case(Isel)
      case(2)
         write(*,*) 'ASSEMBLING MATRIX AND RHS'
      case(1)
         write(*,*) 'ASSEMBLING RHS ONLY'
      case default
         write(*,*) 'ASSEMBLY CASE NOT SUPPORTED!'
         stop
      end select
      ! start_time = MPI_Wtime()
   endif
! MODE_ASSEMBLY = Isel
!
!..assemble global stiffness matrix
!..loop through elements
!
! $OMP PARALLEL                                  &
! $OMP DEFAULT(NONE)                             &
! $OMP SHARED(Isel,RHS,NRELES,ELEM_ORDER,        &
! $OMP        ELEM_NNZ,MAXDOFM,MAXDOFS,NR_PHYSA, &
! $OMP        NR_RHS,ISTC_FLAG,MUMPS_PAR,CLOC   )&
! $OMP PRIVATE(nrdofs,nrdofm,nrdofc,nodm,nrnodm, &
! $OMP         ndofmH,ndofmE,ndofmV,ndofmQ,      &
! $OMP         i,j,k,k1,k2,l,mdle,nod,ndof)
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
   ! if (.not.allocated(LCON0)) 
   ! allocate(LCON0(MAXDOFM))
   allocate(ZLOAD(MAXDOFM))
   allocate(ZTEMP(MAXDOFM**2))
!
! $OMP DO                 &
! $OMP SCHEDULE(DYNAMIC)  &
! $OMP REDUCTION(+:RHS)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      ! write(*,*) 'CLOC(iel)%con=',CLOC(iel)%con
      ! call pause
      if (ISTC_FLAG) then
         call celem_systemI(iel,mdle,2, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      else
         call celem(mdle,2, nrdofs,nrdofm,nrdofc,nodm,  &
            ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,ZLOAD,ZTEMP)
      endif
!
!  ...determine local to global dof connectivities
!       l=0
! !  ...H1 dof
!       do i = nrnodm,1,-1
!          nod = nodm(i)
!          do j=1,ndofmH(i)
!             l=l+1
!             LCON0(l) = NFIRST_DOF(nod)+j
!          enddo
!       enddo
! !  ...H(curl) dof
!       do i = nrnodm,1,-1
!          nod = nodm(i)
!          do j=1,ndofmE(i)
!             l=l+1
!             LCON0(l) = NFIRST_DOF(nod)+ndofmH(i)+j
!          enddo
!       enddo
! !  ...H(div) dof
!       do i = nrnodm,1,-1
!          nod = nodm(i)
!          do j=1,ndofmV(i)
!             l=l+1
!             LCON0(l) = NFIRST_DOF(nod)+ndofmH(i)+ndofmE(i)+j
!          enddo
!       enddo
! !  ...L2 dof
!       if (.not. ISTC_FLAG) then
!          nod = nodm(nrnodm)
!          do j=1,ndofmQ(nrnodm)
!             l=l+1
!             LCON0(l) = NFIRST_DOF(nod)+ndofmH(nrnodm)+ndofmE(nrnodm)+ndofmV(nrnodm)+j
!          enddo
!       endif
!
!  ...number of element (interface) dof
      ! ndof = l

      ndof = CLOC(iel)%ni

      ! write(*,*) 'iel=',iel
      ! write(*,*) 'CLOC'
      ! write(*,*) CLOC(iel)%ni
      ! write(*,*) CLOC(iel)%con(1:ndof)

      ! if (.not.allocated(CLOC(iel)%con)) then
      !    write(*,*)'CLOC(iel)%con is not allocated'
      ! else
      !    write(*,*) 'ubound(CLOC(iel)%con,1)',ubound(CLOC(iel)%con,1)
      ! endif

      ! LCON0(1:ndof) = CLOC(iel)%con(1:ndof)

      ! write(*,*) LCON0
      ! write(*,*) 'ZLOAD'
      ! write(*,*) ZLOAD(1:3)
      ! write(*,*) 'ZTEMP'
      ! write(*,*) ZTEMP(1:3)
      ! write(*,*) ZTEMP(4:6)
      ! write(*,*) ZTEMP(7:9)
!
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      do k1=1,ndof
!     ...global dof is:
         ! i = LCON0(k1)
         i = CLOC(iel)%con(k1)
!     ...Assemble global load vector
         if (isnan(ZLOAD(k1))) then
            write(*,*) 'NaN found in local load vector'
            write(*,*) 'iel,k1=',iel,k1; stop
         endif
         RHS(i) = RHS(i) + ZLOAD(k1)
!     ...Assemble matrix if Isel = 2 only
         if (Isel.ne.2) cycle
!     ...loop through dof `to the right'
         do k2=1,ndof
!        ...global dof is:
            ! j = LCON0(k2)
            j = CLOC(iel)%con(k2)
!        ...assemble
!        ...note: repeated indices are summed automatically by MUMPS
            ! ELEM_NNZ(iel) = ELEM_NNZ(iel) + 1
            k = (k1-1)*ndof + k2
            if (isnan(ZTEMP(k))) then
               write(*,*) 'NaN found in local stiffness matrix'
               write(*,*) 'iel,k1,k2=',iel,k1,k2; stop
            endif
            MUMPS_PAR%A(  ELEM_NNZ(iel)+int8(k)) = ZTEMP(k)
            MUMPS_PAR%IRN(ELEM_NNZ(iel)+int8(k)) = i
            MUMPS_PAR%JCN(ELEM_NNZ(iel)+int8(k)) = j

 !            write(*,*) 'k1,k2,i,j=',k1,k2,i,j
 !            write(*,*) 'ELEM_NNZ(iel)+int8(k)=',ELEM_NNZ(iel)+int8(k)

 !            write(*,*) 'IRN JCN A'
 !            write(*,3993) i , j , ZTEMP(k)
 ! 3993 format(i4,i4,es22.15)

         enddo
      enddo
      ! write(*,*) 'RHS,after iel=',iel
      ! write(*,*) RHS(1:NRDOF_CON)
!
      ! CLOC(iel)%ni = ndof
      ! allocate(CLOC(iel)%con(ndof))
      ! CLOC(iel)%con = LCON0(1:ndof)
!..end of loop through elements
   enddo
! $OMP END DO
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
   deallocate(ZBMOD,ZAMOD,ZLOAD,ZTEMP)
   ! if (allocated(LCON0)) 
   ! deallocate(LCON0)
! $OMP END PARALLEL
!
 !   if (IPRINT_TIME .eq. 1) then
 !      end_time = MPI_Wtime()
 !      Mtime(2) =  end_time-start_time
 !      write(*,1004) Mtime(2)
 ! 1004 format(' STEP 2 finished: ',f12.5,'  seconds',/)
 !   endif
!


end subroutine








subroutine factorize_tangent_mumps
use mumps,     only: MUMPS_PAR,mumps_destroy
!
!
!..MUMPS analysis
   mumps_par%JOB = 1
!
   ! if (IPRINT_TIME .eq. 1) time_stamp = MPI_Wtime()
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
  35 continue
!
!..MUMPS factorization
   mumps_par%JOB = 2
!
   ! if (IPRINT_TIME .eq. 1) time_stamp = MPI_Wtime()
#if C_MODE
   call zmumps(mumps_par)
#else
   call dmumps(mumps_par)
#endif
   if (mumps_par%INFO(1) .ne. 0) then
      write(*,*) 'factorization: mumps_par%INFO(1) .ne. 0'
      if (mumps_par%INFO(1) .eq. -9) then
         write(*,*) 'Increasing workspace, trying factorization again...'
         mumps_par%icntl(14) = mumps_par%icntl(14) + 10 ! increase workspace by 10 percent
         goto 35
      endif
      call mumps_destroy
      stop
   endif
 !   if (IPRINT_TIME .eq. 1) then
 !      time_stamp = MPI_Wtime()-time_stamp
 !      write(*,3002) time_stamp
 ! 3002 format(' - Factorize: ',f12.5,'  seconds')
 !      write(*,1100) '   - memory used in GB    = ',mumps_par%INFO(22)/1000.d0
 !   endif

end subroutine









subroutine get_delta_factor
implicit none
integer :: maxiter_g,iter_g
real*8  :: g0,gs,smin,pm

smin= 0.01d0

if (LINESEARCH_FLAG) then
!!!!!
   write(*,*) 'get_delta_factor: beggining line search...'
   maxiter_g = 5
   iter_g = 1
   LINESEARCH_FACTOR = 0.d0
   call compute_g(0, g0 )
   LINESEARCH_FACTOR = 1.d0
   call compute_g(1, gs )
   do while (abs(gs).gt.abs(g0)*LINESEARCH_TOL.and.iter_g.lt.maxiter_g)
      LINESEARCH_FACTOR =  LINESEARCH_FACTOR * g0 / (g0-gs) 
      call compute_g(1, gs)
      ! gs = gs*LINESEARCH_FACTOR
      iter_g = iter_g + 1
   enddo

   ! if (abs(LINESEARCH_FACTOR).lt.smin) LINESEARCH_FACTOR = smin!sign(smin,LINESEARCH_FACTOR)

   write(*,*) 'get_delta_factor: iter_g, LINESEARCH_FACTOR',iter_g, LINESEARCH_FACTOR
! 
else
   LINESEARCH_FACTOR = 1.d0
endif

end subroutine






! from Matthies and Strang

subroutine get_step_illinois
implicit none
integer :: maxiter_g,iter_g
real*8  :: sa,sb,g0,g1,ga,gb,smax,smin,step, gab

if (LINESEARCH_FLAG) then
!!!!!
   write(*,*) 'get_step_illinois: beggining line search...'
   maxiter_g = 5
   smax= 16.d0
<<<<<<< HEAD
   smin=-1.d0
=======
   ! smin=-16.d0
>>>>>>> c2c6bcb (Meshes, problem setup and run scripts for foam simulation in UW and PR)
   LINESEARCH_FACTOR = 0.d0
   call compute_g(0, g0 )
   LINESEARCH_FACTOR = 1.d0
   call compute_g(1, g1 )

   sa = 1.d0; sb = 0.d0
   ga = g1  ; gb = g0

   do while (ga*gb.gt.0.d0 .and. sa.lt.smax)
      ! gab = ga - gb
      ! if (gab*ga.gt.0.d0) then
      !    sb = min(-1.d0,sb*2.d0)
      !    LINESEARCH_FACTOR = sb; call compute_g(1,gb)
      ! else
      !    sa = sa*2.d0
      !    LINESEARCH_FACTOR = sa; call compute_g(1,ga)
      ! endif
      sb=sa; sa=2.d0*sa; gb=ga
      LINESEARCH_FACTOR = sa; call compute_g(1,ga)
   enddo

   step = sa; g1 = ga

   iter_g = 1
   do while ( ga*gb.lt.0.d0.and.                                &
             (abs(g1).gt.abs(g0)*LINESEARCH_TOL.or.             &
              abs(sb-sa).gt.LINESEARCH_TOL*0.5d0*(sb+sa) ).and. &
             iter_g.lt.maxiter_g                               )

      step = sa - ga*(sa-sb)/(ga-gb)

      LINESEARCH_FACTOR = step; call compute_g(1,g1)

      if (g1*ga.gt.0.d0) then
         gb = 0.5d0*gb
      else
         sb = sa; gb = ga
      endif

      sa = step; ga = g1

      iter_g = iter_g + 1

   enddo

   LINESEARCH_FACTOR = step

   write(*,*) 'get_step_illinois: iter_g, LINESEARCH_FACTOR',iter_g, LINESEARCH_FACTOR
! 
else
   LINESEARCH_FACTOR = 1.d0
endif

end subroutine







subroutine compute_g(Nopt,Gs)
implicit none
integer,intent(in ) :: Nopt
real*8 ,intent(out) :: Gs

! assemble the global residual vector with the current value of LINESEARCH_FACTOR

   select case(Nopt)
   case (0) ! use RES_GLOBAL stored in memory


   case (1) ! obtain a new RES_GLOBAL
      call dealloc_schur
      call assembly_mumps(1,NRDOF_CON,RES_GLOBAL)
      ! write(*,*) 'RES_GLOBAL'
      ! write(*,*) RES_GLOBAL(1:18)
   end select

   Gs = dot_product(RES_GLOBAL,MUMPS_PAR%RHS)

end subroutine







subroutine alloc_BFGS(Neq,Nupdate)
implicit none
integer, intent(in) :: Neq,Nupdate
if (Nupdate.lt.2 .or. Neq.lt.1) return

allocate( V_BFGS(Neq,Nupdate-1),W_BFGS(Neq,Nupdate-1),   &
          RES_PREV(Neq),RES_CURR(Neq)                   )

V_BFGS   = 0.d0
W_BFGS   = 0.d0
RES_PREV = 0.d0
RES_CURR = 0.d0

end subroutine








subroutine dealloc_BFGS

if (allocated(V_BFGS)) deallocate(V_BFGS)

if (allocated(W_BFGS)) deallocate(W_BFGS)

if (allocated(RES_PREV)) deallocate(RES_PREV)

if (allocated(RES_CURR)) deallocate(RES_CURR)

end subroutine








subroutine BFGS_UPDATES(Iter,Nupdate,Neq,Dd,Res0,Res1)
implicit none
integer,intent(in) :: Iter,Nupdate,Neq
real*8 ,intent(in) :: Dd(Neq),Res0(Neq),Res1(Neq)

integer :: m
real*8  :: a,b,c,dres(Neq), tol

tol = 1.d-15

m = mod(Iter-1,Nupdate)+1

if (m.le.1) return

dres = Res1 - Res0

b = dot_product(Dd,dres)
c = dot_product(Res0,Dd)

if (abs(b).gt.tol .and. abs(c).gt.tol) then

   a = max( 0.d0 , - b / c * LINESEARCH_FACTOR )

   V_BFGS(:,m-1) = Dd / b
      
   W_BFGS(:,m-1) = - dres + sqrt(a) * Res0

   write(*,*) 'BFGS_UPDATES: update vectors determined'
else

   V_BFGS(:,m-1) = 0.d0

   W_BFGS(:,m-1) = 0.d0

   write(*,*) 'BFGS_UPDATES: update vectors omitted'

endif



end subroutine







subroutine BFGS_RES(Iter,Nupdate,Neq,Res1,Res_up)
implicit none
integer,intent(in) :: Iter,Nupdate,Neq
real*8 ,intent(in) :: Res1(Neq)
real*8 ,intent(out):: Res_up(Neq)

real*8 :: c
integer:: m,i

m = mod(Iter-1,Nupdate)+1

if (m.le.1) return

Res_up = Res1

do i=m,2,-1
   c      = dot_product(V_BFGS(:,i-1),Res_up)
   Res_up = Res_up + c * W_BFGS(:,i-1)
enddo

end subroutine








subroutine BFGS_SOL(Iter,Nupdate,Neq,Dd)
implicit none
integer,intent(in)    :: Iter,Nupdate,Neq
real*8 ,intent(inout) :: Dd(Neq)

real*8 :: c
integer:: m,i

m = mod(Iter-1,Nupdate)+1

if (m.le.1) return

do i=2,m
   c  = dot_product(W_BFGS(:,i-1),Dd)
   Dd = Dd + c * V_BFGS(:,i-1)
enddo

end subroutine








subroutine solve_delta_mumps(Res)
use mumps,     only: MUMPS_PAR, mumps_destroy
implicit none
VTYPE,intent(in) :: Res(*)

   MUMPS_PAR%RHS(1:NRDOF_CON) = Res(1:NRDOF_CON)
!
!..MUMPS solve
   mumps_par%JOB = 3
!
   ! if (IPRINT_TIME .eq. 1) time_stamp = MPI_Wtime()
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
 !   if (IPRINT_TIME .eq. 1) then
 !      time_stamp = MPI_Wtime()-time_stamp
 !      write(*,3003) time_stamp
 ! 3003 format(' - Solve    : ',f12.5,'  seconds')
 !   endif

end subroutine








subroutine close_solver
   use mumps,       only: mumps_destroy
   use stc,         only: stc_dealloc
   use assembly_sc, only: NFIRST_DOF
   use assembly,    only: MAXDOFS
!
   deallocate(NFIRST_DOF)
   deallocate(MAXDOFS)
   deallocate(ELEM_NNZ)
   call stc_dealloc
!
!..Destroy the instance (deallocate internal data structures)
   call mumps_destroy
!
 !   if (IPRINT_TIME .ge. 1) then
 !      write(*,*)
 !      write(*,1013) sum(Mtime(1:4))
 ! 1013 format(' mumps_sc FINISHED: ',f12.5,'  seconds',/)
 !   endif
!
end subroutine









subroutine resid_norm(Res,Nrdof,Res_norm)
implicit none
real*8, intent(in) :: Res(*)
integer,intent(in) :: Nrdof
real*8, intent(out):: Res_norm

   Res_norm = sqrt(dot_product(Res(1:Nrdof),Res(1:Nrdof)))

end subroutine








subroutine map_dofs_to_local(Vec_global)
   use data_structure3D, only: NRELES
   use assembly,         only: NR_RHS, ZERO
   use stc,              only: CLOC
   use assembly_sc,      only: ZSOL_LOC
implicit none
real*8, intent(in) :: Vec_global(*)

integer :: ndof,iel,i,k1

   ndof = 0
!$OMP PARALLEL
!$OMP DO REDUCTION(MAX:ndof)
   do iel=1,NRELES
      if (CLOC(iel)%ni > ndof) ndof = CLOC(iel)%ni
   enddo
!$OMP END DO
   allocate(ZSOL_LOC(ndof))
!$OMP DO PRIVATE(i,k1) SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      ZSOL_LOC=ZERO
      do k1=1,CLOC(iel)%ni
         i = CLOC(iel)%con(k1)
         ZSOL_LOC(k1) = Vec_global(i)
      enddo
      ! deallocate(CLOC(iel)%con)
      call solout(iel,ndof,NR_RHS,ZERO,ZSOL_LOC)

   enddo
!$OMP END DO
   deallocate(ZSOL_LOC)
!$OMP END PARALLEL

   call dealloc_schur

end subroutine








subroutine dealloc_schur
use data_structure3D, only: NRELES
use stc,              only: CLOC
implicit none
integer :: iel

!$OMP PARALLEL
!$OMP DO
do iel=1,NRELES
   if (allocated(CLOC(Iel)%ASchur)) deallocate(CLOC(Iel)%ASchur)
   if (allocated(CLOC(Iel)%BSchur)) deallocate(CLOC(Iel)%BSchur)
enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine






! subroutine map_new_dofs_to_global(Vec_global)

! !  ...get most recent dofs of approximate solution correctly assembled
!    call assembly_mumps(0,Vec_global)

! end subroutine

! subroutine map_old_dofs_to_global(Vec_global)

! !  ...get previous dofs of approximate solution correctly assembled
!    call assembly_mumps(-1,Vec_global)

! end subroutine

! subroutine get_only_res_global(Vec_global)

! !  ...get rhs (N-R residual) correctly assembled, and not the stiffness matrix
!    call assembly_mumps(1,Vec_global)

! end subroutine
















subroutine shift_com
!
use parameters
use data_structure3D
!
implicit none
integer :: nod,ndofH,ndofE,ndofV,ndofQ

!$OMP DO
   do nod=1,NRNODS
!  ...skip if node inactive
      if (Is_inactive(nod)) cycle
!  ...skip if node out of subdomain
      if (.not.associated(NODES(nod)%dof)) cycle
!  ...retrieve number of dof
      call ndof_nod(NODES(nod)%type,NODES(nod)%order, ndofH,ndofE,ndofV,ndofQ)
! 
!  ...if H1 dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofH)) then
!     ...copy current first copy into second slot
         NODES(nod)%dof%zdofH(NRHVAR+1:2*NRHVAR,1:ndofH) = NODES(nod)%dof%zdofH(1:NRHVAR,1:ndofH)
!     ...reset first copy
         NODES(nod)%dof%zdofH(1:NRHVAR,1:ndofH) = ZERO
      endif
! 
!  ...if H(curl) dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofE)) then
!     ...copy current first copy into second slot
         NODES(nod)%dof%zdofE(NREVAR+1:2*NREVAR,1:ndofE) = NODES(nod)%dof%zdofE(1:NREVAR,1:ndofE)
!     ...reset first copy
         NODES(nod)%dof%zdofE(1:NREVAR,1:ndofE) = ZERO
      endif
! 
!  ...if H(div) dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofV)) then
!     ...copy current first copy into second slot
         NODES(nod)%dof%zdofV(NRVVAR+1:2*NRVVAR,1:ndofV) = NODES(nod)%dof%zdofV(1:NRVVAR,1:ndofV)
!     ...reset first copy
         NODES(nod)%dof%zdofV(1:NRVVAR,1:ndofV) = ZERO
      endif
! 
!  ...if L2 dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofQ)) then
!     ...copy current first copy into second slot
         NODES(nod)%dof%zdofQ(NRQVAR+1:2*NRQVAR,1:ndofQ) = NODES(nod)%dof%zdofQ(1:NRQVAR,1:ndofQ)
!     ...reset first copy
         NODES(nod)%dof%zdofQ(1:NRQVAR,1:ndofQ) = ZERO
      endif
   enddo
!$OMP END DO

end subroutine








subroutine reset_coms
!
use parameters
use data_structure3D
!
implicit none
integer :: nod,ndofH,ndofE,ndofV,ndofQ

!$OMP DO
   do nod=1,NRNODS
!  ...skip if node inactive
      if (Is_inactive(nod)) cycle
!  ...skip if node out of subdomain
      if (.not.associated(NODES(nod)%dof)) cycle
!  ...retrieve number of dof
      call ndof_nod(NODES(nod)%type,NODES(nod)%order, ndofH,ndofE,ndofV,ndofQ)
! 
!  ...if H1 dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofH)) then
!     ...reset both copies
         NODES(nod)%dof%zdofH(1:2*NRHVAR,1:ndofH) = ZERO
      endif
! 
!  ...if H(curl) dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofE)) then
!     ...reset both copies
         NODES(nod)%dof%zdofE(1:2*NREVAR,1:ndofE) = ZERO
      endif
! 
!  ...if H(div) dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofV)) then
!     ...reset both copies
         NODES(nod)%dof%zdofV(1:2*NRVVAR,1:ndofV) = ZERO
      endif
! 
!  ...if L2 dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofQ)) then
!     ...reset both copies
         NODES(nod)%dof%zdofQ(1:2*NRQVAR,1:ndofQ) = ZERO
      endif
   enddo
!$OMP END DO

end subroutine








subroutine update_sol(Factor)
!
use parameters
use data_structure3D
!
implicit none
real*8, intent(in) :: Factor

integer :: nod,ndofH,ndofE,ndofV,ndofQ

!$OMP DO PRIVATE( ndofH,ndofE,ndofV,ndofQ)
   do nod=1,NRNODS
!  ...skip if node inactive
      if (Is_inactive(nod)) cycle
!  ...skip if node out of subdomain
      if (.not.associated(NODES(nod)%dof)) cycle
!  ...retrieve number of dof
      call ndof_nod(NODES(nod)%type,NODES(nod)%order, ndofH,ndofE,ndofV,ndofQ)
! 
!  ...if H1 dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofH)) then
         ! write(*,*) 'nod',nod
         ! write(*,*) 'before'
         ! write(*,*) NODES(nod)%dof%zdofH(:,1)
!     ...copy current first copy into second slot
         NODES(nod)%dof%zdofH(NRHVAR+1:2*NRHVAR,1:ndofH) = NODES(nod)%dof%zdofH(NRHVAR+1:2*NRHVAR,1:ndofH) &
                                                         + Factor * NODES(nod)%dof%zdofH(1:NRHVAR,1:ndofH)
!     ...reset first copy
         NODES(nod)%dof%zdofH(1:NRHVAR,1:ndofH) = ZERO
         ! write(*,*) 'nod',nod
         ! write(*,*) 'after'
         ! write(*,*) NODES(nod)%dof%zdofH(:,1)
      endif
! 
!  ...if H(curl) dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofE)) then
!     ...copy current first copy into second slot
         NODES(nod)%dof%zdofE(NREVAR+1:2*NREVAR,1:ndofE) = NODES(nod)%dof%zdofE(NREVAR+1:2*NREVAR,1:ndofE) &
                                                         + Factor * NODES(nod)%dof%zdofE(1:NREVAR,1:ndofE)
!     ...reset first copy
         NODES(nod)%dof%zdofE(1:NREVAR,1:ndofE) = ZERO
      endif
! 
!  ...if H(div) dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofV)) then
!     ...copy current first copy into second slot
         NODES(nod)%dof%zdofV(NRVVAR+1:2*NRVVAR,1:ndofV) = NODES(nod)%dof%zdofV(NRVVAR+1:2*NRVVAR,1:ndofV) &
                                                         + Factor * NODES(nod)%dof%zdofV(1:NRVVAR,1:ndofV)
!     ...reset first copy
         NODES(nod)%dof%zdofV(1:NRVVAR,1:ndofV) = ZERO
      endif
! 
!  ...if L2 dof are allocated, proceed with shifting
      if (associated(NODES(nod)%dof%zdofQ)) then
!     ...copy current first copy into second slot
         NODES(nod)%dof%zdofQ(NRQVAR+1:2*NRQVAR,1:ndofQ) = NODES(nod)%dof%zdofQ(NRQVAR+1:2*NRQVAR,1:ndofQ) &
                                                         + Factor * NODES(nod)%dof%zdofQ(1:NRQVAR,1:ndofQ)
!     ...reset first copy
         NODES(nod)%dof%zdofQ(1:NRQVAR,1:ndofQ) = ZERO
      endif
   enddo
!$OMP END DO

end subroutine

!------------------------------------------------------------------------------------
!
!  routine: copy_coms
!
!  last modified: Jan 2019
!
!  purpose: copies all solution dofs from one component set to another
!
!  input:   - No1: component to copy from
!           - No2: component to copy to
!
!------------------------------------------------------------------------------------
!
subroutine copy_coms(No1,No2)
!
   use parameters
   use data_structure3D
!
   implicit none
!
   integer, intent(in)  :: No1,No2
!
   integer :: nod, nf, nt, nn2, i
!
!------------------------------------------------------------------------------------
!
!..check consistency
   if ((No1.lt.0).or.(No2.lt.0).or.(No1.gt.NRCOMS).or.(No2.gt.NRCOMS)) then
      write(*,*) 'copy_coms: No1,No2,NRCOMS = ', No1,No2,NRCOMS, ' . stop.'
      stop
   endif
   if (No1.eq.No2) return
!
!..loop through active nodes
!$OMP PARALLEL DO          &
!$OMP PRIVATE(nf,nt,nn2,i) &
!$OMP SCHEDULE(DYNAMIC)
   do nod=1,NRNODS
      if (Is_inactive(nod)) cycle
      if (.not. associated(NODES(nod)%dof)) cycle
!
!  ...H1 dof
      if (.not. associated(NODES(nod)%dof%zdofH)) goto 10
      nf = (No1-1)*NRHVAR
      nt = (No2-1)*NRHVAR
      nn2 = ubound(NODES(nod)%dof%zdofH,2)
      if(nn2.gt.0) then
         do i=1,NRHVAR
            NODES(nod)%dof%zdofH(nt+i,1:nn2) = NODES(nod)%dof%zdofH(nf+i,1:nn2)
         enddo
      endif
  10  continue
!
!  ...H(curl) dof
      if (.not. associated(NODES(nod)%dof%zdofE)) goto 20
      nf = (No1-1)*NREVAR
      nt = (No2-1)*NREVAR
      nn2 = ubound(NODES(nod)%dof%zdofE,2)
      if(nn2.gt.0) then
         do i=1,NREVAR
            NODES(nod)%dof%zdofE(nt+i,1:nn2) = NODES(nod)%dof%zdofE(nf+i,1:nn2)
         enddo
      endif
  20  continue
!
!  ...H(div) dof
      if (.not. associated(NODES(nod)%dof%zdofV)) goto 30
      nf = (No1-1)*NRVVAR
      nt = (No2-1)*NRVVAR
      nn2 = ubound(NODES(nod)%dof%zdofV,2)
      if(nn2.gt.0) then
         do i=1,NRVVAR
            NODES(nod)%dof%zdofV(nt+i,1:nn2) = NODES(nod)%dof%zdofV(nf+i,1:nn2)
         enddo
      endif
  30  continue
!
!  ...L2 dof
      if (.not. associated(NODES(nod)%dof%zdofQ)) goto 40
      nf = (No1-1)*NRQVAR
      nt = (No2-1)*NRQVAR
      nn2 = ubound(NODES(nod)%dof%zdofQ,2)
      if(nn2.gt.0) then
         do i=1,NRQVAR
            NODES(nod)%dof%zdofQ(nt+i,1:nn2) = NODES(nod)%dof%zdofQ(nf+i,1:nn2)
         enddo
      endif
  40  continue
!
!..end of loop through nodes
   enddo
!$OMP END PARALLEL DO
!
end subroutine copy_coms
!
!
end module