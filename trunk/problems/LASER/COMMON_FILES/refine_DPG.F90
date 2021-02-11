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
subroutine refine_DPG(Irefine,Nreflag,Factor,Nflag,PhysNick,Ires, Nstop)
!
   use commonParam
   use control
   use data_structure3D
   use environment  , only: QUIET_MODE
   use assembly_sc  , only: NRDOF_CON, NRDOF_TOT
   use parametersDPG, only: NORD_ADD
   use par_mesh     , only: DISTRIBUTED,HOST_MESH
   use mpi_param    , only: ROOT,RANK,NUM_PROCS
   use MPI          , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD,   &
                            MPI_REAL8,MPI_INTEGER,MPI_IN_PLACE,MPI_MAX
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real(8), intent(in)  :: Factor
   integer, intent(in)  :: Nflag(NR_PHYSA)
   integer, intent(in)  :: PhysNick
   logical, intent(in)  :: Ires
   integer, intent(out) :: Nstop
!
   integer, parameter :: max_step = 20
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
!
!..adaptive refinements
   integer :: nr_elem_ref,ref_pml
   integer :: ref_ndom(4)
   integer :: mdle_ref(NRELES)
   integer :: mdle_pml(NRELES)
   real(8) :: res_sum
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
!
   real(8) :: resid_subd,resid_tot,elem_resid_max
   real(8) :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
   real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd
!
   integer :: i,ic,mdle,iel,kref,subd,count,ierr
   real(8) :: x(3), xnod(3,8)
!
!..element type
   character(len=4) :: etype
!
   real(8) :: MPI_Wtime,start_time,end_time
!
!..printing flag
   integer :: iprint = 0
   !character(len=8) :: filename
!
!-----------------------------------------------------------------------
!
!..initialize
   elem_resid(1:NRELES)    = 0.d0
   elem_ref_flag(1:NRELES) = 0
   mdle_ref(1:NRELES)      = 0
   mdle_pml(1:NRELES)      = 0
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
      if (USE_PML .and. Is_pml(mdle)) then
         ! treat PML differently if needed
         mdle_pml(iel) = 1
      else
         ! outside PML
      endif
      if(ires) then
         call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
         resid_subd = resid_subd + elem_resid(iel)
         if (elem_resid(iel) > elem_resid_max) elem_resid_max = elem_resid(iel)
      endif
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
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
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
   if (USE_PML) then
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         if (mdle_pml(iel) .gt. 0) then
            NODES(mdle)%visit = 1
         else
            NODES(mdle)%visit = 0
         endif
      enddo
   endif
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (.not. QUIET_MODE) then
      if (RANK .eq. ROOT) write(*,1010) end_time-start_time
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
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
!
!..mark elements for adaptive refinement
   if (Irefine .eq. IADAPTIVE) then
      nr_elem_ref = 0
      ! strategy 1 (greedy strategy)
!      do iel=1,NRELES
!         mdle = ELEM_ORDER(iel)
!         call find_domain(mdle, i)
!         if(elem_resid(iel)>0.75d0*elem_resid_max) then
!            nr_elem_ref = nr_elem_ref + 1
!            iel_ref(nr_elem_ref) = iel
!            !if (RANK.eq.ROOT) write(*,1020) 'ndom=',i,', mdle=',mdle,', r =',elem_resid(iel)
!         endif
 1020    format(A,I2,A,I9,A,es12.5)
!      enddo
      ! strategy 2 (Doerfler's strategy)
      ! sort elements (iel_ref) according to residual values
      mdle_ref(1:NRELES) = ELEM_ORDER(1:NRELES)
      call qsort_duplet(mdle_ref,elem_resid,NRELES,1,NRELES)
      ! verify the array is sorted
      do iel=1,NRELES-1
         if (elem_resid(iel) .lt. elem_resid(iel+1)) then
            write(*,1021) 'elem_resid not sorted: elem_resid(',iel  ,') = ',elem_resid(iel), &
                                               ', elem_resid(',iel+1,') = ',elem_resid(iel+1)
 1021       format(A,I9,A,es12.5,A,I9,A,es12.5,/)
            stop
         endif
      enddo
      !
      if (RANK.eq.ROOT) write(*,*) 'Array is sorted.'
      if (RANK.eq.ROOT) write(*,*) 'The following elements are marked for refinement:'
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
      !
      res_sum = 0.d0
      ref_ndom(1:4)=0; ref_pml=0 ! stats
      do iel=1,NRELES
         mdle = mdle_ref(iel)
         call find_domain(mdle, i)
         ref_ndom(i) = ref_ndom(i)+1
         !
         if (USE_PML .and. NODES(mdle)%visit .eq. 1) then
            ! treat PML differently if needed
            ref_pml = ref_pml+1
            !if (RANK.eq.ROOT) write(*,1022) 'PML: ndom=',i,', mdle=',mdle,', res =',elem_resid(iel)
         else
            ! outside PML
            !if (RANK.eq.ROOT) write(*,1022) '     ndom=',i,', mdle=',mdle,', res =',elem_resid(iel)
         endif
 1022    format(A,I2,A,I9,A,es12.5)
         nr_elem_ref = nr_elem_ref + 1
         res_sum = res_sum + elem_resid(iel)
         if(res_sum > 0.5d0*resid_tot) exit
      enddo
      ! WRITE HISTORY FILE FOR DEBUGGING
!      filename='HIST.dat'
!      open(UNIT=9,FILE=filename,FORM="FORMATTED",ACCESS="APPEND",STATUS="UNKNOWN",ACTION="WRITE")
!      write(UNIT=9, FMT="(I6)") nr_elem_ref
!      do iel=1,nr_elem_ref
!         write(UNIT=9, FMT="(I6)") mdle_ref(iel)
!      enddo
!      close(UNIT=9)
      ! END WRITE FOR DEBUGGING
      if (RANK.eq.ROOT) write(*,1023) 'nr_elem_ref = ', nr_elem_ref
      if (RANK.eq.ROOT) write(*,1023) '        PML = ', ref_pml
      if (RANK.eq.ROOT) write(*,1023) '     ndom 1 = ', ref_ndom(1)
      if (RANK.eq.ROOT) write(*,1023) '     ndom 2 = ', ref_ndom(2)
      if (RANK.eq.ROOT) write(*,1023) '     ndom 3 = ', ref_ndom(3)
      if (RANK.eq.ROOT) write(*,1023) '     ndom 4 = ', ref_ndom(4)
 1023 format(A,I7)
   endif
!
!  END MARKING OF ELEMENTS FOR ADAPTIVE REFINEMENTS
!
   if (Irefine.eq.ICORE .or. Irefine.eq.ICLAD) then
      nr_elem_ref = 0
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         call find_domain(mdle, i)
         if (Irefine.eq.ICORE .and. (i.eq.1 .or. i.eq.2)) then
            nr_elem_ref = nr_elem_ref + 1
            mdle_ref(nr_elem_ref) = mdle
         !elseif (Irefine.eq.ICLAD .and. (i.eq.3 .or. i.eq.4)) then
         elseif (Irefine.eq.ICLAD .and. (i.eq.3)) then
            nr_elem_ref = nr_elem_ref + 1
            mdle_ref(nr_elem_ref) = mdle
         endif
      enddo
   endif
!
!..use appropriate strategy depending on refinement type
   Nstop = 0
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
!  ...anisotropic refinements
      case(IANISOTROPIC)
!     ...anisotropic href (z)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call global_href_aniso(0,1) ! refine in z
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call update_gdof
         call update_Ddof
!  ...adaptive refinements
      case(IADAPTIVE,ICORE,ICLAD)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         if (RANK .eq. ROOT)  write(*,*) 'Starting adaptive refinements...'
         if (RANK .eq. ROOT)  write(*,*) 'NRELES = ', NRELES
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
            etype = NODES(mdle)%type
            select case(etype)
               case('mdlb')
                  !kref = 111 ! iso
                  !kref = 110  ! radial
                  kref = 10 ! refining in r
                  !kref = 100 ! refining in theta
               case('mdlp')
                  !kref = 11  ! iso
                  kref = 10   ! radial
               case default
                  write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
                  call pause
            end select
            call refine(mdle,kref)
         enddo
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call close_mesh
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time
!         call enforce_min_rule
!         call enforce_max_rule ! MAX RULE ISSUE NEEDS TO BE FIXED
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
  2020 format(' refinement : ',f12.5,'  seconds')
  2025 format(' close mesh : ',f12.5,'  seconds')
  2030 format(A,I8,', ',I9)
!
!-----------------------------------------------------------------------
!                              FINALIZE
!-----------------------------------------------------------------------
!
   if (IBCFLAG .eq. 3) then
      call propagate_flag(2,9)
      call propagate_flag(3,9)
   endif
!
   90 continue
!
end subroutine refine_DPG
!
!
!-----------------------------------------------------------------------
! subroutine: href_solve (uniform refinements and solve)
!-----------------------------------------------------------------------
subroutine href_solve(Nflag,PhysNick, Nstop)
!
   use commonParam
   use control
   use data_structure3D
   use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
   use mpi_param     , only: RANK,ROOT
   use par_mesh      , only: DISTRIBUTED,HOST_MESH,distr_mesh
   use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_partition
!
   implicit none
!
   integer, intent(in)  :: Nflag(NR_PHYSA)
   integer, intent(in)  :: PhysNick
   integer, intent(out) :: Nstop
!
   integer :: nsteps,count,src,ierr,i
   logical :: ires
!
   Nstop = 0
   ires = .true.
!
   if (RANK .eq. ROOT) then
      nsteps=0
      do while (nsteps.le.0)
         write(*,*) 'Provide: number of uniform h-refinements'
         read(*,*) nsteps
      enddo
   else
      write(6,405) '[', RANK, '] : ','Waiting for broadcast from master...'
 405  format(A,I4,A,A)
   endif
   count = 1; src = ROOT
   call MPI_BCAST (nsteps,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
!
!..initial solve
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      call par_mumps_sc('H')
   else
      !call mumps_sc('G')
      call pardiso_sc('H')
   endif
!
!..do refinements and solve
   do i=0,nsteps
!  ...display error and refine if necessary
      if (i.ne.nsteps) then
!     ...Uniform refinement and solve
         call refine_DPG(IUNIFORM,1,0.25d0,Nflag,PhysNick,ires, Nstop)
         if (DISTRIBUTED .and. (.not. HOST_MESH)) then
            !call zoltan_w_set_lb(1)
            !call distr_mesh
            !call print_partition
            call par_mumps_sc('H')
         else
            !call mumps_sc('H')
            call pardiso_sc('H')
         endif
      else
!     ...Last step only display (no refinement)
         call refine_DPG(INOREFINEMENT,1,0.25d0,Nflag,PhysNick,ires, Nstop)
      endif
   enddo
!
end subroutine href_solve
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
