!-----------------------------------------------------------------------
! routine name - Hp adapt
!---------------------------------------------------------------------

!last revision: 18 Nov 2022

!Purpose: Main driver for Hp adaptation: full fledged

! arguments:





!---------------------------------------------------------

subroutine HpAdapt

    use common_prob_data
    use control
    use data_structure3D
    use environment     , only: QUIET_MODE
    use assembly_sc     , only: NRDOF_CON, NRDOF_TOT
    use parametersDPG
    use par_mesh        , only: DISTRIBUTED, HOST_MESH
    use mpi_param       , only: ROOT, RANK,NUM_PROCS
    use MPI             , only: MPI_COMM_WORLD, MPI_SUM,MPI_COMM_WORLD, &
                                MPI_REAL8, MPI_INTEGER, MPI_IN_PLACE, MPI_MAX


    implicit none
    integer, parameter :: adap_strat = 1 !0 is for greedy strat and 1 is for Doerfler 
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
    real(8) :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
    real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd


    integer :: i,ic,mdle,iel,kref,subd,count,ierr
    real(8) :: x(3), xnod(3,8)
    integer :: nord_new,is,nord,nordx,nordy,nordz,naux,pord
    !for the time being for Hex only
    integer :: nr_sons,kref_out
    integer :: first_son
    integer :: mdle_sons(8)
    real(8) :: error_org,error_p,error_hopt
    real(8), allocatable :: xnod_ref(:,:,:)

!..element type
    character(len=4) :: etype
!
    real(8) :: MPI_Wtime,start_time,end_time
!
!..printing flag
    integer :: iprint = 0

! to hold the index of sons of hp refined elements on coarse mesh to obtain fine mesh.

!..initialize
    elem_resid(1:NRELES)    = 0.d0
    elem_ref_flag(1:NRELES) = 0
    mdle_ref(1:NRELES)      = 0


    Nflag(1:NR_PHYSA) = (/0,0,1,1/) 


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
   call MPI_ALLREDUCE(error_subd,error_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(rnorm_subd,rnorm_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
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
7005  format(' mesh |','     NE     |','  nrdof_tot |','  nrdof_con |','    residual   |','   residual rate  ',/)
7006  format(' mesh |','     NE     |','  nrdof_tot |','  nrdof_con |','    residual   |','   residual rate  |', &
           ' field error  |','rel field error|','   error rate ',/)
!
do i=1,istep
   if (NEXACT.eq.0) then
      write(*,7003) i,nelem_mesh(i),nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i)
7003    format(2x,i2,'  | ',3(i10,' | '),es12.5,'  |',f7.2)
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
             if(elem_resid(iel) > Factor * elem_resid_max) then
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
             if(res_sum > Factor*resid_tot) exit
 
          enddo
          if (RANK.eq.ROOT) write(*,1023) 'nr_elem_ref = ', nr_elem_ref
          1023 format(A,I7)
       endif
endif

allocate(xnod_ref(nr_elem_ref,3,MAXbrickH))
!..use appropriate strategy depending on refinement type
Nstop = 0

! select case(Irefine)

!     case(IADAPTIVE)
!         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!         if (RANK .eq. ROOT)  write(*,*) 'Starting hp refining the mesh to compute fine mesh solution'
!         if (RANK .eq. ROOT)  write(*,*) 'NRELES = ', NRELES        

!         if(adap_strat .eq. 0) then
!             do iel = 1,nr_elem_ref
!              mdle = mdle_ref(iel)
!              etype = NODES(mdle)%type
!              nord = NODES(mdle)%order
!              select case(etype)
!              case('mdlb')
!                 kref = 111
!                 nord_new = nord + 111
                
!                 call decode(nord_new,naux,nordz)
!                 call decode(naux,nordx,nordy)
!                 pord = MAX(nordx,nordy,nordz)
!                 if (pord .gt. MAXP) then
!                     write(*,1001) 'local pref: mdle,p,MAXP = ',mdle,pord,MAXP,'. stop.'
!                     stop
!                     1001 format(A,I7,I3,I3,A)
!                 endif

!              case default
!                 write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
!                 call pause
!              end select
!             !p refine and then h-refine
!              call nodmod(mdle,nord_new)
!              call refine(mdle,kref)
            
!             enddo
!         elseif(adap_strat .eq. 1) then
!             do iel = 1,nr_elem_ref
!                 mdle = mdle_ref(iel)
!                 etype = NODES(mdle)%type
!                 nord = NODES(mdle)%order
!                 select case(etype)
!                    case('mdlb')
!                       kref = 111 ! iso
!                       nord_new = nord + 111
                
!                       call decode(nord_new,naux,nordz)
!                       call decode(naux,nordx,nordy)
!                       pord = MAX(nordx,nordy,nordz)
!                      !  write(*,*) mdle,pord,MAXP
!                       if (pord .gt. MAXP) then
!                           write(*,1002) 'local pref: mdle,p,MAXP = ',mdle,pord,MAXP,'. stop.'
!                           stop
!                           1002 format(A,I7,I3,I3,A)
!                       endif   
!                      !  call nodmod(mdle,nord_new)                  
!                       !kref = 110  ! radial
!                       ! kref = 10 ! refining in r
!                       !kref = 100 ! refining in theta
!                    ! case('mdlp')
!                    !    !kref = 11  ! iso
!                    !    kref = 10   ! radial
!                    case default
!                       write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
!                       call pause
!                 end select
!             !p refine and then h-refine
!                !  call nodcor(mdle,dummy_xnod)
!                !  xnod_ref(iel,1:3,1:MAXbrickH) = dummy_xnod(1:3,1:MAXbrickH)
                
!                 call nodmod(mdle,nord_new)
!                !  call refine(mdle,kref)
!                !  call break(mdle,kref)
!             enddo
!         endif
!       !   call enforce_min_rule
!         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
!         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!         call close_mesh
!         call enforce_min_rule
!         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time
!   !..raise order of approximation on non-middle nodes by enforcing minimum rule
        
!         call par_verify
!         call update_gdof
!         call update_Ddof

!         if (RANK .eq. ROOT) write(*,*) 'NRELES = ', NRELES
!         if (RANK .eq. ROOT) write(*,*) 'Finished hp refining selected elements'


!     case default; Nstop = 1
! end select
! 2020 format(' fine mesh refinement : ',f12.5,'  seconds')
! 2025 format(' close mesh : ',f12.5,'  seconds')
! 2030 format(A,I8,', ',I9)

call Finegrid_padap(nr_elem_ref,mdle_ref,xnod_ref)

! solving on the hp refined mesh to fine mesh solution

if (DISTRIBUTED .and. (.not. HOST_MESH)) then
   !call zoltan_w_set_lb(1)
   !call distr_mesh
   !call print_partition
   call par_mumps_sc('G')
else
   call mumps_sc('G')
endif

! write(*,*) NRDOF_TOT
call exact_error
! fine mesh solver and now looping over elements to make choices
do iel = 1,nr_elem_ref  
   mdle = mdle_ref(iel)
   etype = NODES(mdle)%type
   ! nord = NODES(mdle)%order
   ! nr_sons = NODES(mdle)%nr_sons
   ! first_son = NODES(mdle)%first_son
   ! write(*,*) "number of sons = ",nr_sons,mdle
   select case(etype)

      case('mdlb')
         !storing mdle for element sons of elements refined on the coarse mesh
         ! to get the fine mesh
         ! do is = 1,8
         !    mdle_sons(is) = first_son + is - 1
         ! enddo
         call project_p(mdle,error_org,error_p)
         call project_h(Mdle,error_org,error_hopt,kref_out)    
         ! call project_h(mdle,error_h_best)
      case default
         write(*,*) "Element type not recongnized"

   end select



enddo


end subroutine HpAdapt