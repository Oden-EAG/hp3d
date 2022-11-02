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
subroutine refine_DPG

    use common_prob_data
    use control
    use data_structure3D
    use environment  , only: QUIET_MODE
    use assembly_sc  , only: NRDOF_CON, NRDOF_TOT
    use parametersDPG, only: NORD_ADD
    use par_mesh     , only: DISTRIBUTED,HOST_MESH
    use mpi_param    , only: ROOT,RANK,NUM_PROCS
    use MPI          , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD,   &
                             MPI_REAL8,MPI_INTEGER,MPI_IN_PLACE,MPI_MAX
    
   !  integer, intent(in)  :: Irefine
   !  integer, intent(in)  :: Nreflag
   !  real(8), intent(in)  :: Factor
   !  integer, intent(in)  :: Nflag(NR_PHYSA)
   !  integer, intent(in)  :: PhysNick
   !  logical, intent(in)  :: Ires
   !  integer, intent(out) :: Nstop
    
    integer, parameter :: adap_strat = 0 !0 is for greedy strat and 1 is for doefeler
    ! parameters below were input to the function before but currently used as finxed parameters for debug
    integer, parameter :: physNick = 1 !if exact then adapting to reduce in error L2 solution u
    integer, parameter :: max_step = 20
    integer, parameter :: Factor = 0.75
    logical :: Ires = .true.
    integer, parameter :: Irefine = 2
    
    
    integer :: Nflag(NR_PHYSA)
    integer :: Nreflag
    integer :: Nstop

    integer, dimension(max_step), save :: nrdof_tot_mesh
    integer, dimension(max_step), save :: nrdof_con_mesh
    integer, dimension(max_step), save :: nelem_mesh
    real(8), dimension(max_step), save :: residual_mesh
    real(8), dimension(max_step), save :: rate_mesh
    real(8), dimension(max_step), save :: error_mesh
    real(8), dimension(max_step), save :: rel_error_mesh
    real(8), dimension(max_step), save :: rate_error_mesh
 !
   real(8) :: elem_resid(NRELES)
   integer :: elem_ref_flag(NRELES)


!..adaptive refinements
   integer :: nr_elem_ref
   integer :: ref_ndom(4)
   integer :: mdle_ref(NRELES)
   real(8) :: res_sum
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
   real(8) :: resid_subd,resid_tot,elem_resid_max
   real(8) :: errorQ,rnormQ 
   real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd


   integer :: i,ic,mdle,iel,kref,subd,count,ierr
   real(8) :: x(3), xnod(3,8)

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


   Nflag(1:NR_PHYSA) = (/0,0,1,0/) 

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


!..initialize residual and error
   resid_subd = 0.d0
   error_subd = 0.d0
   rnorm_subd = 0.d0
   elem_resid_max = 0.d0   

!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()

   !..residual/error computation
   !
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(errorH,errorE,errorV,errorQ,            &
!$OMP         rnormH,rnormE,rnormV,rnormQ,            &
!$OMP         mdle,subd)                              &
!$OMP SCHEDULE(DYNAMIC)                               &
!$OMP REDUCTION(+:resid_subd,error_subd,rnorm_subd)   &
!$OMP REDUCTION(MAX:elem_resid_max)
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      if (DISTRIBUTED .and. (RANK .ne. subd)) cycle

      if(ires) then ! compute the max residual and the residual in the subdomain
         call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
         resid_subd = resid_subd + elem_resid(iel) 
         if(elem_resid(iel) > elem_resid_max) elem_resid_max = elem_resid(iel)
      endif

      if(NEXACT .ge. 1) then
         call element_error(mdle,Nflag,                     &
                         errorH,errorE,errorV,errorQ,    &
                         rnormH,rnormE,rnormV,rnormQ)
         select case(PhysNick)
            case(1)
            error_subd = error_subd + errorQ
            rnorm_subd = rnorm_subd + rnormQ
         end select
      endif                   

   enddo
!$OMP END PARALLEL DO

resid_tot = 0.d0; error_tot = 0.d0; rnorm_tot = 0.d0
if(DISTRIBUTED) then !MPI reduction op to ollect the error and residual
   count = 1
   call MPI_ALLREDUCE(resid_subd,resid_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(error_subd,resid_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(rnorm_subd,resid_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid_max,count,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr) !max over all elements
   count = NRELES
   call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   ! call MPI_ALLREDUCE(MPI_IN_PLACE,mdle_pml,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
else
   resid_tot = resid_subd
   error_tot = error_subd
   rnorm_tot = rnorm_subd
endif

!..end timer
call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
if (.not. QUIET_MODE) then
   if (RANK .eq. ROOT) write(*,1010) end_time-start_time
1010 format(' elem_residual  : ',f12.5,'  seconds')
endif

if (RANK .ne. ROOT) goto 70

!..update and display convergence history
nelem_mesh(istep) =  NRELES
residual_mesh(istep) = sqrt(resid_tot)
if (NEXACT.ge.1) then
   error_mesh(istep)     = sqrt(error_tot)
   rel_error_mesh(istep) = sqrt(error_tot/rnorm_tot)
endif

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

!..print out the history of refinements
write(*,*)
write(*,*) 'HISTORY OF REFINEMENTS'
if (NEXACT.eq.0) write(*,7005)
if (NEXACT.gt.0) write(*,7006)
7005  format(' mesh |','  nrdof_tot |','  nrdof_con |','    residual   |','   residual rate  ',/)
7006  format(' mesh |','     NE     |','  nrdof_tot |','  nrdof_con |','    residual   |','   residual rate  |', &
           ' field error  |','rel field error|','   error rate ',/)
!
do i=1,istep
   if (NEXACT.eq.0) then
      write(*,7003) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i)
7003    format(2x,i2,'  | ',2(i10,' | '),es12.5,'  |',f7.2)
   else
      write(*,7004) i,nelem_mesh(i),nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i), &
                    error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
7004    format(2x,i2,'  | ',3(i10,' | '),es12.5,'  |',f7.2,'          ', &
             2(' | ',es12.5),'  |',f7.2)
   endif
   if (i .eq. istep) write(*,*)
enddo

70 continue

!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
!..mark elements for adaptive refinement
if (Irefine .eq. IADAPTIVE) then
   nr_elem_ref = 0
   ! strategy 1 (greedy strategy)
      if(adap_strat .eq. 0) then
         do iel = 1,NRELES
            mdle = ELEM_ORDER(iel)
            call find_domain(mdle,i)
            if(elem_resid(iel) > 0.8d0 * elem_resid_max) then
               nr_elem_ref = nr_elem_ref + 1
               mdle_ref(nr_elem_ref) = mdle
               if (RANK.eq.ROOT) write(*,1020) 'ndom=',i,', mdle=',mdle,', r =',elem_resid(iel)
            endif
            1020    format(A,I2,A,I9,A,es12.5)
         enddo
      elseif(adap_strat .eq. 1) then    ! strategy 2: Doerfler's Strategy

         mdle_ref(1:NRELES) = ELEM_ORDER(1:NRELES)
         call qsort_duplet(mdle_ref,elem_resid,NRELES,1,NRELES)
         ! verify the array is sorted
         do iel=1,NRELES-1
            if (elem_resid(iel) .lt. elem_resid(iel+1)) then
               write(*,1021) 'elem_resid not sorted: elem_resid(',iel  ,') = ',elem_resid(iel), &
                                                ', elem_resid(',iel+1,') = ',elem_resid(iel+1)
            1021     format(A,I9,A,es12.5,A,I9,A,es12.5,/)
               stop
            endif
         enddo
         !
         if (RANK.eq.ROOT) write(*,*) 'Array is sorted.'
         if (RANK.eq.ROOT) write(*,*) 'The following elements are marked for refinement:'
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)      

         res_sum = 0.d0
         ref_ndom(1:4)=0

         do iel = 1,NRELES
            mdle = mdle_ref(iel)
            call find_domain(mdle,i)
            ref_ndom(i) = ref_ndom(i) + 1

            nr_elem_ref = nr_elem_ref + 1
            res_sum = res_sum + elem_resid(iel)
            if(res_sum > 0.5d0*resid_tot) exit

         enddo
         if (RANK.eq.ROOT) write(*,1023) 'nr_elem_ref = ', nr_elem_ref
         1023 format(A,I7)
      endif
endif

!  END MARKING OF ELEMENTS FOR ADAPTIVE REFINEMENTS

!..use appropriate strategy depending on refinement type
Nstop = 0

select case(Irefine)
   case(IUNIFORM)
!     isotropic href
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call global_href
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
      if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
      call update_gdof
      call update_Ddof


   case(IADAPTIVE)
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      if (RANK .eq. ROOT)  write(*,*) 'Starting adaptive refinements...'
      if (RANK .eq. ROOT)  write(*,*) 'NRELES = ', NRELES

      if(adap_strat .eq. 0) then
        do iel = 1,nr_elem_ref
         mdle = mdle_ref(iel)
         etype = NODES(mdle)%type
         select case(etype)
         case('mdlb')
            kref = 111
         case default
            write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
            call pause
         end select
         call refine(mdle,kref)
        enddo
      elseif(adap_strat .eq. 1) then
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
            etype = NODES(mdle)%type
            select case(etype)
               case('mdlb')
                  kref = 111 ! iso
                  !kref = 110  ! radial
                  ! kref = 10 ! refining in r
                  !kref = 100 ! refining in theta
               ! case('mdlp')
               !    !kref = 11  ! iso
               !    kref = 10   ! radial
               case default
                  write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
                  call pause
            end select
            call refine(mdle,kref)
         enddo
      endif
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
      if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call close_mesh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
      if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time

      call par_verify
      call update_gdof
      call update_Ddof
!
      if (RANK .eq. ROOT) write(*,*) 'NRELES = ', NRELES
      if (RANK .eq. ROOT) write(*,*) 'Finished adaptive refinements...'

   case default; Nstop = 1
end select
2020 format(' refinement : ',f12.5,'  seconds')
2025 format(' close mesh : ',f12.5,'  seconds')
2030 format(A,I8,', ',I9)
!
!-----------------------------------------------------------------------
!                              FINALIZE
!-----------------------------------------------------------------------
!
if (IBCFLAG .eq. 3) then
   call propagate_flag(3,3)
   call propagate_flag(5,3)
endif
!
90 continue


end subroutine refine_DPG


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



subroutine adap_solve()
   !
      use common_prob_data
      use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
      use mpi_param     , only: RANK,ROOT
      use par_mesh      , only: DISTRIBUTED,HOST_MESH,distr_mesh
      use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_partition
   !
      implicit none
   !
      integer :: nsteps,count,src,ierr,i
   !
      if (RANK .eq. ROOT) then
         nsteps=0
         do while (nsteps.le.0)
            write(*,*) 'Provide: number of adaptive h-refinements'
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
         call par_mumps_sc('G')
      else
         call mumps_sc('G')
      endif
   !
   !..do refinements and solve
      do i=0,nsteps
   !  ...display error and refine if necessary
         if (i.ne.nsteps) then
   !     ...adaptive refinement and solve
            call refine_DPG
            if (DISTRIBUTED .and. (.not. HOST_MESH)) then
               !call zoltan_w_set_lb(1)
               !call distr_mesh
               !call print_partition
               call par_mumps_sc('G')
            else
               call mumps_sc('G')
            endif
         else
   !     ...Last step only display (no refinement)
            call refine_DPG
         endif
      enddo
   !
   end subroutine adap_solve
   