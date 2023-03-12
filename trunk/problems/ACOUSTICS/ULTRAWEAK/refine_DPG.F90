!--------------------------------------------------------------------
!
!     routine name      - refine_DPG
!
!--------------------------------------------------------------------
!
!     last revision:    - Oct 2019
!
!     purpose:          - refines elements assuming problem has
!                         already been solved. If uniform refinements
!                         it refines everything. Otherwise it follows
!                         greedy strategy based on residual to
!                         determine which elements to refine.
!                         It also displays dof, residuals and
!                         relevant errors at each step.
!
!     arguments:
!
!     in:
!             Irefine   = 0 if no refinement -> only display
!                       = 1 if uniform refinements
!                       = 2 if anisotropic refinements (provide kref)
!                       = 3 if adaptive refinements
!             Nreflag   = 1 for h-refinements
!                       = 2 for p-refinements
!                       = 3 for h- or p-refinement, depending
!                           upon the element size
!             Factor    - if element error \ge Factor*max_error
!                         then the element is refined
!             Nflag     - integer vector of length NR_PHYSA to indicate
!                         which PHYS Attribute one should compute the error for
!             PhysNick  - nickname for PHYS Attribute: physNick = HEVQ
!                         e.g., physNick = 1001 means 1H1 and 1L2 variable, etc.
!             Ires      - logical, .true. if compute residual is true
!     out:
!             Nstop     = 1 if refinement did not increase dof
!
!---------------------------------------------------------------------
!
subroutine refine_DPG(Irefine,Nreflag,Factor, Nstop)
!
   use common_prob_data_UW
   use control
   use element_data
   use data_structure3D
   use environment  , only: QUIET_MODE
   use assembly_sc  , only: NRDOF_CON, NRDOF_TOT
   use parametersDPG, only: NORD_ADD
   use par_mesh     , only: DISTRIBUTED,HOST_MESH
   use mpi_param    , only: ROOT,RANK,NUM_PROCS
   use MPI          , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD,   &
                            MPI_REAL8,MPI_INTEGER,MPI_IN_PLACE,MPI_MAX
   use mg_balance   , only: UNIFORM_REF, MDLE_REF, nr_elem_ref,      &
                            ref_type, ref_kind
   use mg_stats     , only: GSTAT
   use mg_data_structure, only: CURRENT_GRID
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real(8), intent(in)  :: Factor
   integer, intent(out) :: Nstop
!
   integer :: Nflag(NR_PHYSA)
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
   integer, allocatable :: elem_pref(:), elem_pref_ord(:)
!
!..adaptive refinements
   real(8) :: res_sum
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
   integer, save :: STAGE = 1
!
   real(8) :: resid_subd,resid_tot,elem_resid_max
   real(8) :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
   real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd
!
   integer :: i,ic,mdle,iel,kref,subd,count,ierr,icpref
   integer :: iref,nord_comp(3),nordx_new,nordy_new,nordz_new
   integer :: nord,nord_new
   integer :: idec_pref,refine_strategy
   real(8) :: x(3), xnod(3,8)
!
!..element type
   integer :: ntype
!
   real(8) :: MPI_Wtime,start_time,end_time
!
!..printing flag
   integer :: iprint = 1
   !character(len=8) :: filename
!
!-----------------------------------------------------------------------
!
!..Set refinement strategy
   refine_strategy = 2
!
!..For MGCG solver we first mark, then exit and refine
   if (STAGE .eq. 2) goto 80
!
!..field variables flag - used inside ndof_DPG and for error calc.
   nflag(1:NR_PHYSA)=(/0,0,1/)
!
!..initialize
   elem_resid(1:NRELES)    = 0.d0
   elem_ref_flag(1:NRELES) = 0
   allocate(MDLE_REF(NRELES)); MDLE_REF = 0
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
      if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
!
      call elem_residual_fi(mdle, elem_resid(iel),elem_ref_flag(iel))
!
      resid_subd = resid_subd + elem_resid(iel)
      if (elem_resid(iel) > elem_resid_max) elem_resid_max = elem_resid(iel)
      if (NEXACT.ge.1) then
         call element_error(mdle,Nflag, errorH,errorE,errorV,errorQ,  &
                                        rnormH,rnormE,rnormV,rnormQ)
         error_subd = error_subd + errorH + errorE + errorV + errorQ
         rnorm_subd = rnorm_subd + rnormH + rnormE + rnormV + rnormQ
      endif
   enddo
!$OMP END PARALLEL DO
!
   resid_tot = 0.d0; error_tot = 0.d0; rnorm_tot = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_ALLREDUCE(resid_subd,resid_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(error_subd,error_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rnorm_subd,rnorm_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid_max,count,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      count = NRELES
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   else
      resid_tot = resid_subd
      error_tot = error_subd
      rnorm_tot = rnorm_subd
   endif
!
   GSTAT(CURRENT_GRID)%error = sqrt(error_tot)
   GSTAT(CURRENT_GRID)%resid = sqrt(resid_tot)
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (.not. QUIET_MODE) then
      if (RANK .eq. ROOT) write(*,1010) end_time-start_time
 1010 format('elem_residual  : ',f12.5,'  seconds')
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
 7005  format(' mesh |',' nrdof_tot |',' nrdof_con |','    residual   |','   residual rate  ',/)
 7006  format(' mesh |',' nrdof_tot |',' nrdof_con |','    residual   |','   residual rate  |', &
              ' field error  |','rel field error|','   error rate ',/)
!
   do i=1,istep
      if (NEXACT.eq.0) then
         write(*,7003) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i)
 7003    format(2x,i2,'  | ',2(i9,' | '),es12.5,'  |',f7.2)
      else
         write(*,7004) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i), &
                       error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
 7004    format(2x,i2,'  | ',2(i9,' | '),es12.5,'  |',f7.2,'          ', &
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
            if(elem_resid(iel)>0.5d0*elem_resid_max) then
               nr_elem_ref = nr_elem_ref + 1
               MDLE_REF(nr_elem_ref) = mdle
               !if (RANK.eq.ROOT) write(*,1020) 'ndom=',i,', mdle=',mdle,', r =',elem_resid(iel)
            endif
 1020       format(A,I2,A,I9,A,es12.5)
         enddo
      else
!
!     Strategy 2 (Doerfler's strategy)
!
!     ...sort elements (iel_ref) according to residual values
         MDLE_REF(1:NRELES) = ELEM_ORDER(1:NRELES)
         call qsort_duplet(mdle_ref,elem_resid,NRELES,1,NRELES)
!     ...verify the array is sorted
         do iel=1,NRELES-1
            if (elem_resid(iel) .lt. elem_resid(iel+1)) then
               write(*,1021) char(9), 'elem_resid not sorted: elem_resid(',iel  ,') = ',elem_resid(iel), &
                                         ', elem_resid(',iel+1,') = ',elem_resid(iel+1)
1021           format(A,I9,A,es12.5,A,I9,A,es12.5,/)
               stop
            endif
         enddo
!
         if (RANK.eq.ROOT) write(*,*) char(9), 'Array is sorted.'
         if (RANK.eq.ROOT) write(*,*) char(9), 'The following elements are marked for refinement:'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
         res_sum = 0.d0
         do iel=1,NRELES
            mdle = MDLE_REF(iel)
!
1022        format(A,I2,A,I9,A,es12.5)
            nr_elem_ref = nr_elem_ref + 1
            res_sum = res_sum + elem_resid(iel)
            if(res_sum > 0.75d0*resid_tot) exit
         enddo
!
         if (RANK.eq.ROOT) write(*,1023) char(9), 'nr_elem_ref = ', nr_elem_ref
1023     format(A,A,I7)
      endif
   endif
!
   Nstop = 0
!..decide refinement types
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
         UNIFORM_REF = .true.
!
!  ...adaptive refinements
      case(IADAPTIVE)
         UNIFORM_REF = .false.
!
!     ...allocate refinement arrays
         allocate(ref_type(nr_elem_ref)); ref_type(:) = 0
         allocate(ref_kind(nr_elem_ref)); ref_kind(:) = 0
!
!     ...get refinement kind
         do iel = 1,nr_elem_ref
            mdle = MDLE_REF(iel)
!
!        ...Select_refinements can only be called on SUBD
!           Nreflag=1,2 doesn't need to cycle, but communication is cheap anyway
            call get_subd(mdle, subd)
            if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
!
            ntype = NODES(mdle)%ntype
!
!        ...set potential h-ref
            select case(ntype)
               case(MDLB)
                  kref = 111 ! iso
               case(MDLP)
                  kref = 11  ! iso
               case default
                  write(*,*) char(9), 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',ntype
                  stop 1
            end select
!
!        ...set potential p-ref
            nord = NODES(mdle)%order
            select case(ntype)
            case(MDLB)
               call decod(nord,10,3,nord_comp(1:3))
               nordx_new = min(nord_comp(1)+1,MAXP)
               nordy_new = min(nord_comp(2)+1,MAXP)
               nordz_new = min(nord_comp(3)+1,MAXP)
               nord_new = nordx_new*100+nordy_new*10+nordz_new
            case default
               write(*,*) 'refine_DPG: UNFINISHED' ; stop 1
            end select

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
               ref_type(iel) = 2
               ref_kind(iel) = nord_new
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
                     ref_type(iel) = 0
                     ref_kind(iel) = 0
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
   if (STAGE .eq. 1) then
!  ...set stage to refine so that next call will execute refinements
      STAGE = 2
      goto 90
   endif
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
   80 continue
!
   Nstop = 0
!..use appropriate strategy depending on refinement type
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
!     ...isotropic href
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call global_href
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call update_gdof
         call update_Ddof
!  ...adaptive refinements
      case(IADAPTIVE)
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
            mdle = MDLE_REF(iel)
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
               if (RANK.eq.ROOT) write(*,*) 'Elem:', mdle, 'not refined, ref_type = ',ref_type(iel)
            end select
         enddo
!
!     ...deallocate refinement type arrays
         deallocate(ref_type,ref_kind)
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call close_mesh
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time


         select case(Nreflag)
         case(2)
            write(*,*) 'refine_DPG: ENFORCING MAX RULE --- UNFINISHED'
!            call enforce_max_rule
            stop 1
         case(3)
            write(*,*) 'refine_DPG: ENFORCING MAX RULE --- UNFINISHED'
!            call enforce_max_rule
            stop 1
         case(4)
            if (icpref.gt.0) then
               call execute_pref(elem_pref(1:icpref), icpref)
               call enforce_min_rule
            endif
         end select
!
         deallocate(elem_pref,elem_pref_ord)
!
!         call enforce_min_rule
!         call enforce_max_rule ! MAX RULE ISSUE NEEDS TO BE FIXED
         call par_verify
         call update_gdof
         call update_Ddof
!
         if (RANK .eq. ROOT) write(*,*) char(9), 'NRELES = ', NRELES
         if (RANK .eq. ROOT) write(*,*) char(9), 'Finished adaptive refinements...'
!
      case default; Nstop = 1
!
   end select
!
!..Reset STAGE for next refinement cycle
   deallocate(MDLE_REF)
   STAGE = 1
!
  2020 format(' refinement : ',f12.5,'  seconds')
  2025 format(' close mesh : ',f12.5,'  seconds')
  2030 format(A,I8,', ',I9)
!
   if (IBC_PROB.eq.3 .or. IBC_PROB.eq.4 .or. IBC_PROB.eq.6) call propagate_flag(2,3)
!
   90 continue
!
end subroutine refine_DPG
!
!
!-----------------------------------------------------------------------
!
!    routine name:      - qsort_duplet
!
!-----------------------------------------------------------------------
!
!    latest revision:   - Oct 2019
!
!    purpose:           - sorts an array of duplets (iel,residual) with
!                         residual (sort key) in descending order
!                         (initial call needs: First = 1, Last = N)
!
!    arguments:
!           in/out
!                       - Iel_array: 1D integer array (element indices)
!                       - Residuals: 1D real    array (residual values)
!           in
!                       - N        : size of array
!                       - First    : first index of current partition
!                       - Last     : last  index of current partition
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
recursive subroutine qsort_duplet(Iel_array,Residuals,N,First,Last)
!
   implicit none
!..declare variables
   integer , intent(in)    :: N,First,Last
   integer , intent(inout) :: Iel_array(N)
   real(8) , intent(inout) :: Residuals(N)
!
   real(8) :: x,pivot
   integer :: i,j,l
!
   pivot = Residuals((First+Last) / 2)
   i = First
   j = Last
!..iterate through the array to be sorted
   do
!  ...find first element from the left that needs to be swapped
      do while ((Residuals(i) > pivot))
         i = i + 1
      end do
!  ...find first element from the right that needs to be swapped
      do while ((pivot > Residuals(j)))
         j = j - 1
      end do
!  ...end loop if no elements need to be swapped
      if (i >= j) exit
!  ...swap the elements
      l = Iel_array(i); Iel_array(i) = Iel_array(j); Iel_array(j) = l
      x = Residuals(i); Residuals(i) = Residuals(j); Residuals(j) = x
      i = i + 1
      j = j - 1
   end do
   if (First < i-1) call qsort_duplet(Iel_array,Residuals,N,First,i-1 )
   if (j+1 < Last)  call qsort_duplet(Iel_array,Residuals,N,j+1,  Last)
!
end subroutine qsort_duplet
