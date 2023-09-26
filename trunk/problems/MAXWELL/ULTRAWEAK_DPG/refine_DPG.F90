
!------------------------------------------------------------------------------
!> @brief      Refines mesh based on current solution and provided flags
!!
!> @param[in]  Irefine     - refinement approach
!!                             0: no refinement -> only display
!!                             1: uniform refinements
!!                             2: adaptive refinements
!> @param[in]  Nreflag     - refinement kind
!!                             1: h-refinements
!!                             2: p-refinements
!!                             3: h- or p-refinement
!> @param[in]  Factor      - threshold used in marking refinement
!!
!> @param[out] Nstop       - 1 if refinement did not increase dof; 0 otherwise
!!
!> @date       July 2023
!------------------------------------------------------------------------------
subroutine refine_DPG(Irefine,Nreflag,Factor, Nstop)
!
   use commonParam
   use control
   use data_structure3D
   use environment  , only: QUIET_MODE
   use assembly_sc  , only: NRDOF_CON, NRDOF_TOT
   use parametersDPG, only: NORD_ADD
   use par_mesh     , only: DISTRIBUTED,HOST_MESH
   use mpi_param    , only: ROOT,RANK,NUM_PROCS
   use sorts        , only: qsort_duplet
   use MPI          , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD,   &
                            MPI_REAL8,MPI_INTEGER,MPI_IN_PLACE,MPI_MAX
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real(8), intent(in)  :: Factor
   integer, intent(out) :: Nstop
!
   integer :: Nflag(NR_PHYSA)
   integer :: PhysNick
!
   integer, parameter :: max_step = 50
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real(8), dimension(max_step), save :: residual_mesh
   real(8), dimension(max_step), save :: rate_mesh
   real(8), dimension(max_step), save :: error_mesh
   real(8), dimension(max_step), save :: rel_error_mesh
   real(8), dimension(max_step), save :: rate_error_mesh
!
   real(8) :: elem_resid(NRELES)
   integer :: elem_ref_flag(NRELES)
   integer :: nr_elem_ref
   integer, allocatable :: elem_pref(:), elem_pref_ord(:)
   integer, allocatable :: mdle_ref(:), ref_type(:), ref_kind(:)
!
!..adaptive refinements
   integer :: ref_pml
   integer :: ref_ndom(4)
   integer :: mdle_pml(NRELES)
   real(8) :: res_sum
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
   integer, save :: STAGE = 1
!
   integer :: iref, nord_comp(3), nordx_new, nordy_new, nordz_new
   integer :: nord, nord_new, icpref, idec_pref
!
   real(8) :: errorH, errorE, errorV, errorQ
   real(8) :: rnormH, rnormE, rnormV, rnormQ
!
   real(8) :: resid_subd, resid_tot, elem_resid_max
   real(8) :: error_tot,  rnorm_tot
   real(8) :: error_subd, rnorm_subd
!
   integer :: nod, nel, mdle, iel, kref, subd
   integer :: count, ierr, i, ic
   real(8) :: x(3), xnod(3,8)
   integer :: refine_strategy
!
!..element type
   integer :: ntype
!
   real(8) :: MPI_Wtime,start_time,end_time
!
!..printing flag
   integer :: iprint = 0
!
!-----------------------------------------------------------------------
!
!..initialize
   elem_resid(1:NRELES)    = 0.d0
   elem_ref_flag(1:NRELES) = 0
   mdle_pml(1:NRELES)      = 0
   allocate(mdle_ref(NRELES)); mdle_ref = 0
!
   physNick = 1
   Nflag=0; Nflag(2)=1
!
!..increase step if necessary
!..irefineold=0 means no refinement was made in the previous step
!..if first call or if a refinement was made, increase step
   if ((istep.eq.0).or.(irefineold.ne.0)) istep=istep+1
   irefineold=Irefine
!
!..get dof count from solver
   nrdof_tot_mesh(istep) = NRDOF_TOT
   nrdof_con_mesh(istep) = NRDOF_CON
   if (NRDOF_TOT .eq. 0) then
      istep=istep-1
      goto 70
   endif
!
!..initialize residual and error
   resid_subd = 0.d0
   error_subd = 0.d0
   rnorm_subd = 0.d0
   elem_resid_max = 0.d0
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..residual/error computation
!
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(errorH,errorE,errorV,errorQ,            &
!$OMP         rnormH,rnormE,rnormV,rnormQ,            &
!$OMP         mdle,subd)                              &
!$OMP SCHEDULE(DYNAMIC)                               &
!$OMP REDUCTION(+:resid_subd,error_subd,rnorm_subd)   &
!$OMP REDUCTION(MAX:elem_resid_max)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
!
      if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
!
      call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
!
      resid_subd = resid_subd + elem_resid(iel)
      if (elem_resid(iel) > elem_resid_max) elem_resid_max = elem_resid(iel)
!
      if (NEXACT.ge.1) then
         call element_error(mdle,Nflag, errorH,errorE,errorV,errorQ,  &
                                        rnormH,rnormE,rnormV,rnormQ)
         select case(PhysNick)
            case(1)
               error_subd = error_subd + errorQ
               rnorm_subd = rnorm_subd + rnormQ
            case(10)
               error_subd = error_subd + errorV
               rnorm_subd = rnorm_subd + rnormV
            case(100)
               error_subd = error_subd + errorE
               rnorm_subd = rnorm_subd + rnormE
            case(1000)
               error_subd = error_subd + errorH
               rnorm_subd = rnorm_subd + rnormH
            case(1001)
               error_subd = error_subd + errorH + errorQ
               rnorm_subd = rnorm_subd + rnormH + rnormQ
            case default
               error_subd = error_subd + errorH + errorE + errorV + errorQ
               rnorm_subd = rnorm_subd + rnormH + rnormE + rnormV + rnormQ
         end select
      endif
   enddo
!$OMP END PARALLEL DO
!
   resid_tot = 0.d0; error_tot = 0.d0; rnorm_tot = 0.d0
   if (DISTRIBUTED) then
      count = 1
      call MPI_ALLREDUCE(resid_subd,resid_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(error_subd,error_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rnorm_subd,rnorm_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid_max,count,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      count = NRELES
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,mdle_pml,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
   else
      resid_tot = resid_subd
      error_tot = error_subd
      rnorm_tot = rnorm_subd
   endif
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (.not. QUIET_MODE) then
      if (RANK .eq. ROOT) write(*,1010) end_time - start_time
 1010 format(' elem_residual  : ',f12.5,'  seconds')
   endif
!
   if (RANK .ne. ROOT) goto 70
!
!..update and display convergence history
   residual_mesh(istep) = sqrt(resid_tot)
   if (NEXACT.ge.1) then
      error_mesh(istep)     = sqrt(error_tot)
      rel_error_mesh(istep) = sqrt(error_tot/rnorm_tot)
   endif
!
!..compute decrease rate for the residual and error
   select case(istep)
   case(1)
      rate_mesh(istep) = 0.d0
      if (NEXACT.ge.1) rate_error_mesh(istep) = 0.d0
   case default
      rate_mesh(istep) =  &
      log(residual_mesh(istep-1)/residual_mesh(istep))/  &
      log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      if (NEXACT.ge.1) then
         rate_error_mesh(istep) = &
         log(error_mesh(istep-1)/error_mesh(istep))/  &
         log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      endif
   end select
!
!..print out the history of refinements
   write(*,*)
   write(*,*) 'HISTORY OF REFINEMENTS'
   if (NEXACT.eq.0) write(*,7005)
   if (NEXACT.gt.0) write(*,7006)
 7005  format(' mesh |','  nrdof_tot |','  nrdof_con |','    residual   |','   residual rate  ',/)
 7006  format(' mesh |','  nrdof_tot |','  nrdof_con |','    residual   |','   residual rate  |', &
              ' field error  |','rel field error|','   error rate ',/)
!
   do i=1,istep
      if (NEXACT.eq.0) then
         write(*,7003) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i)
 7003    format(2x,i2,'  | ',2(i10,' | '),es12.5,'  |',f7.2)
      else
         write(*,7004) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i), &
                       error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
 7004    format(2x,i2,'  | ',2(i10,' | '),es12.5,'  |',f7.2,'          ', &
                2(' | ',es12.5),'  |',f7.2)
      endif
      if (i .eq. istep) write(*,*)
   enddo
!
   70 continue
!
!..mark elements for adaptive refinement
   if (Irefine .eq. IADAPTIVE) then
      nr_elem_ref = 0
      if (refine_strategy .eq. 1) then
!
!     Strategy 1 (greedy strategy)
!
         do iel=1,NRELES
            mdle = ELEM_ORDER(iel)
            call find_domain(mdle, i)
            if(elem_resid(iel) > Factor * elem_resid_max) then
               nr_elem_ref = nr_elem_ref + 1
               mdle_ref(nr_elem_ref) = mdle
            endif
 1020       format(A,I2,A,I9,A,es12.5)
        enddo
      else
!
!     Strategy 2 (Doerfler's strategy)
!
!     ...sort elements (iel_ref) according to residual values
         mdle_ref(1:NRELES) = ELEM_ORDER(1:NRELES)
         call qsort_duplet(mdle_ref,elem_resid,NRELES,1,NRELES)
!
!     ...verify the array is sorted
         do iel=1,NRELES-1
            if (elem_resid(iel) .lt. elem_resid(iel+1)) then
               write(*,1021) 'elem_resid not sorted: elem_resid(',iel  ,') = ',elem_resid(iel), &
                                                  ', elem_resid(',iel+1,') = ',elem_resid(iel+1)
 1021       format(A,I9,A,es12.5,A,I9,A,es12.5,/)
               stop
            endif
         enddo
!
         res_sum = 0.d0
         do iel=1,NRELES
            mdle = ELEM_ORDER(iel)
!
            nr_elem_ref = nr_elem_ref + 1
            res_sum = res_sum + elem_resid(iel)
            if(res_sum > Factor * resid_tot) exit
         enddo
      endif
!
      if (RANK.eq.ROOT) write(*,1023) 'nr_elem_ref = ', nr_elem_ref
 1023 format(A,I7)
   endif
!
!..use appropriate strategy depending on refinement type
   Nstop = 0
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
         continue
!
!  ...adaptive refinements
      case(IADAPTIVE)
!
!     ...allocate refinement arrays
         allocate (ref_type(nr_elem_ref)); ref_type(:) = 0
         allocate (ref_kind(nr_elem_ref)); ref_kind(:) = 0
!
!     ...select refinement kind
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
!
            call get_subd(mdle, subd)
            if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
!
            ntype = NODES(mdle)%ntype
!
!        ...set potential h-ref
            select case(ntype)
               case(MDLB)
                  kref = 111     ! iso
                  !kref = 110    ! radial
                  !kref = 10     ! refining in r
                  !kref = 100    ! refining in theta
               case(MDLP)
                  kref = 11      ! iso
                  !kref = 10     ! radial
               case(MDLN,MDLD)
                  call get_isoref(mdle, kref)
               case default
                  write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ', s_type(ntype)
                  call pause
            end select
!
!        ...set potential p-ref
            nord = NODES(mdle)%order
            select case(ntype)
            case(MDLB)
               call decod(nord,10,3,nord_comp(1:3))
               nordx_new = min(nord_comp(1) + 1, MAXP)
               nordy_new = min(nord_comp(2) + 1, MAXP)
               nordz_new = min(nord_comp(3) + 1, MAXP)
               nord_new  = nordx_new*100 + nordy_new*10 + nordz_new
            case(MDLP)
               call decod(nord,10,2,nord_comp(1:2))
               nordx_new = min(nord_comp(1) + 1, MAXP)
               nordy_new = min(nord_comp(2) + 1, MAXP)
               nord_new  = nordx_new*10 + nordy_new
            case(MDLN,MDLD)
               nord_new  = min(nord + 1, MAXP)
            case default
               write(*,*) 'refine_DPG: UNFINISHED' ; stop 1
            end select
!
            if (nord_new.eq.nord) then
               idec_pref=0
            else
               idec_pref=1
            endif
!
            select case(Nreflag)
            case(1)
               ref_type(iel) = 1
               ref_kind(iel) = kref
            case(2)
               if (idec_pref.eq.1) then
                  ref_type(iel) = 2
                  ref_kind(iel) = nord_new
!           ...p-refinement selected but already at max order
               else
                  ref_type(iel) = 0
                  ref_kind(iel) = 0
               endif
            case(3,4)
               call select_refinement(mdle,iref)
               select case(iref)
!           ...h-refinement selected
               case(1)
                  ref_type(iel) = 1
                  ref_kind(iel) = kref
!           ...p-refinement selected
               case(2)
                  if (idec_pref.eq.1) then
                     ref_type(iel) = 2
                     ref_kind(iel) = nord_new
!              ...p-refinement selected but already at max order
                  else
                     ref_type(iel) = 1
                     ref_kind(iel) = kref
                  endif
               end select
            end select
         enddo
!
!     ...communicate refinements
         if (DISTRIBUTED .and. nr_elem_ref .gt. 0) then
            count = nr_elem_ref
            call MPI_ALLREDUCE(MPI_IN_PLACE,ref_type,count,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,ref_kind,count,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
         endif
!
      case default; Nstop = 1
!
   end select
!
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
!..use appropriate strategy depending on refinement type
   Nstop = 0
   select case(Irefine)
!
!  ...uniform refinements
      case(IUNIFORM)
!
!     ...isotropic href
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
         select case (Nreflag)
         case(1)
            call global_href
         case(2)
            call global_pref
         case default
            call global_href
         end select
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call update_gdof
         call update_Ddof
!
!  ...adaptive refinements
      case(IADAPTIVE)
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         if (RANK .eq. ROOT)  write(*,*) char(9), 'Starting adaptive refinements...'
         if (RANK .eq. ROOT)  write(*,*) char(9), 'NRELES = ', NRELES
!
!     ...arrays to buffer p-refinements
         allocate (elem_pref(nr_elem_ref))
         allocate (elem_pref_ord(nr_elem_ref))
         icpref = 0
!
!     ...execute h-refinements, buffer prefinements
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
!
            select case(ref_type(iel))
            case(0)
               cycle
            case(1)
               kref = ref_kind(iel)
               call refine(mdle,kref)
            case(2)
               icpref = icpref + 1
               elem_pref(icpref) = mdle
!           ...not currently used, execute_pref does uniform p ref
               elem_pref_ord(icpref) = nord_new
            case default
               if (RANK.eq.ROOT) write(*,*) 'Elem:', mdle, 'not refined, ref_type = ',ref_type(iel), 'invalid'
            end select
         enddo
!
!     ...deallocate refinement type arrays
         deallocate(ref_type, ref_kind)
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call close_mesh
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time
!
         select case(Nreflag)
         case(2)
            call execute_pref(elem_pref(1:icpref), icpref)
         case(3)
            if (icpref.gt.0) then
               call execute_pref(elem_pref(1:icpref), icpref)
               call enforce_min_rule
            endif
         end select
!
         deallocate(elem_pref,elem_pref_ord)
!
         call par_verify
         call update_gdof
         call update_Ddof
!
         if (RANK .eq. ROOT) write(*,*) 'NRELES = ', NRELES
         if (RANK .eq. ROOT) write(*,*) 'Finished adaptive refinements...'
!
      case default; Nstop = 1
!
   end select
!
   deallocate(mdle_ref)
!
  2020 format(' refinement : ',f12.5,'  seconds')
  2025 format(' close mesh : ',f12.5,'  seconds')
  2030 format(A,I8,', ',I9)
!
!-----------------------------------------------------------------------
!                              FINALIZE
!-----------------------------------------------------------------------
!
   if (IBCFLAG .eq. 3) call propagate_flag(3,3)
!
end subroutine refine_DPG






!------------------------------------------------------------------------------
!> @brief      criterion set by user for choosing between h and p refinement
!!
!> @param[in]  Mdle     - element (middle node) number
!> @param[out] Iref     - refinement kind
!!                         1: h-refinement
!!                         2: p-refinement
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine select_refinement(Mdle,Iref)
!
      use data_structure3D
      use commonParam
      use element_data
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Iref
!
      integer :: i,j,nod,nodp,iv,ie,nrv,nre,nel
      integer :: nodesl(27), norientl(27)
      real*8  :: h, lambda
!
!---------------------------------------------------------------------
!
! TODO: define critereon here
      Iref = 1
!
   end subroutine select_refinement

