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
    use par_mesh        , only: DISTRIBUTED, HOST_MESH, distr_mesh
    use mpi_param       , only: ROOT, RANK,NUM_PROCS
    use MPI             , only: MPI_COMM_WORLD, MPI_SUM,MPI_COMM_WORLD, &
                                MPI_REAL8, MPI_INTEGER, MPI_IN_PLACE, MPI_MAX


    implicit none
    integer, parameter :: adap_strat = 1 !0 is for greedy strat and 1 is for Doerfler 
    ! parameters below were input to the function before but currently used as finxed parameters for debug
    integer, parameter :: physNick = 1 !if exact then adapting to reduce in error L2 solution u
    integer, parameter :: max_step = 40
    real(8), parameter :: Factor = 0.75
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
   !  integer,   allocatable :: Ref_indicator_flags(:,:)
    integer,   allocatable :: Ref_indicator_flags(:)
    real(8),   allocatable :: Nref_grate(:)
    real(8) :: grate_mesh
    integer,   allocatable :: info_old_nodes(:,:)
    integer :: active_node_count
    integer :: NRNODS_coarse,NRNODS_fine
    logical, allocatable   :: act_flags_coarse(:)
    integer :: count_deact_ref
    integer, allocatable   :: deact_ref_index(:)
    integer :: count_add(4)
    integer :: Poly_flag
    integer,   allocatable :: flag_pref(:)
    integer :: hx,hy,hz, nr_sons_intent,nr_sons_close,mdle_child
    integer, allocatable :: pref_intent(:), pref_close(:)  
 !

    integer, save :: istep = 0
    integer, save :: irefineold = 0
    real(8) :: resid_subd,resid_tot,elem_resid_max
    real(8) :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
    real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd
    integer :: NrdofH,NrdofE,NrdofV,NrdofQ


    integer :: i,ic,mdle,iel,kref,subd,count,ierr,inod
    real(8) :: x(3), xnod(3,MAXbrickH)
    integer :: nord_new,is,nord,nordx,nordy,nordz,naux,pord
    !for the time being for Hex only
    integer :: nr_sons,kref_out
    integer :: first_son
    integer :: mdle_sons(8)
    real(8) :: error_org,error_p,error_hopt,rate_p
    real(8), allocatable :: xnod_ref(:,:,:)

!..element type
    character(len=4)  :: etype
    character(len=15) :: mesh_filename
    character(len=4) :: rankno

! for coarse mesh retrieval
    type(node), allocatable :: NODES_cp(:)
    integer,    allocatable :: ELEM_ORDER_cp(:)
    integer,    allocatable :: ELEM_SUBD_cp(:)
    integer :: NRELES_old,NRNODS_old,NRDOFSH_old,NRDOFSE_old,NRDOFSV_old,NRDOFSQ_old
    integer :: NPNODS_old,MAXNODS_old,NRELIS_old
    integer :: kref_close, kref_intent,href_count
    integer,   allocatable :: write_kref(:,:),write_pref(:,:)
!
    real(8) :: MPI_Wtime,start_time,end_time, timer_a, timer_b,timer_c, net_time_a, net_time_b
    real(8), dimension(max_step,2), save :: timer_save
!
!..printing flag
    integer :: iprint = 0

! to hold the index of sons of hp refined elements on coarse mesh to obtain fine mesh.

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

   !  write(*,*) "istep = ",istep

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
        call get_subd(mdle, subd) ! get the owner process or subdomain of the element.
        if (DISTRIBUTED .and. (RANK .ne. subd)) cycle !dont compute if element doesnt belong to the process.

        if(ires) then ! compute the max residual and the residual in the subdomain
           call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel)) !returns the square of the residual
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
if(DISTRIBUTED) then !MPI reduction op to collect the error and residual
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

open(22,file="error.dat",status='replace')
! write(4,*) nr_elem_ref
   do i  = 1,istep

      write(22,*) i,nelem_mesh(i),nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i), &
                 error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)

   enddo
close(22)

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
          if (RANK.eq.ROOT) write(*,*) 'nr_elem_ref = ', nr_elem_ref,res_sum,Factor*resid_tot
         !  1023 format(A,I7)
       endif
endif

allocate(xnod_ref(nr_elem_ref,3,MAXbrickH))
!..use appropriate strategy depending on refinement type
Nstop = 0

call MPI_BARRIER (MPI_COMM_WORLD, ierr)

!start copying the current nodes structure into a copy nodes structure
allocate(NODES_cp(MAXNODS))
call Nodes_copy(NODES_cp)
NRELES_old = NRELES
NRNODS_old = NRNODS
NRDOFSH_old = NRDOFSH
NRDOFSE_old = NRDOFSE
NRDOFSV_old = NRDOFSV
NRDOFSQ_old = NRDOFSQ
NPNODS_old = NPNODS
MAXNODS_old = MAXNODS
NRELIS_old = NRELIS

!----------------------------------------------------------------------------
! call dumpout_hp3d("befor_mesh_copy")

call MPI_BARRIER (MPI_COMM_WORLD, ierr)

allocate(flag_pref(nr_elem_ref)) !flag_pref contains the flags if element has been hp refined
flag_pref = ZERO
!hp refining the marked elements
call Finegrid_padap(nr_elem_ref,mdle_ref,xnod_ref,flag_pref)


! solving on the hp refined mesh to fine mesh solution

if (DISTRIBUTED .and. (.not. HOST_MESH)) then
   !call zoltan_w_set_lb(1)
   !call distr_mesh
   !call print_partition
   call par_mumps_sc('G')
else
   call mumps_sc('G')
endif
call exact_error

call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!--------------------------------------------------------------------------
! if(NUM_PROCS .lt. 10) then
   allocate(Ref_indicator_flags(nr_elem_ref * 10))
   Ref_indicator_flags = ZERO
   !first column will contain indicator flag for
   ! kind of refinement: 1 for h-ref and 2 for p-ref, 2nd column will have the refinement flag: eg: 000 means no
   ! h-ref and only p-ref and non-zero href indicator means href, from 3rd column onwards we have the polynomial
   ! distribution for each subson generated during h refinement, Note: the memory allocation for p-dist of subson
   ! is done for isotropic brick element as it will produce maximum subsons and thus, memory alocation will work for 
   !every other element like prism, pyramid and tets.

   allocate(Nref_grate(nr_elem_ref))
   Nref_grate = ZERO
   net_time_a = 0.0
   net_time_b = 0.0

   !$OMP PARALLEL DO                                          &
   !$OMP PRIVATE(mdle,etype,subd,error_org,rate_p,poly_flag) &
   !$OMP SCHEDULE(DYNAMIC)
      do iel = 1,nr_elem_ref  
         mdle = mdle_ref(iel)
         etype = NODES(mdle)%type        
         call get_subd(mdle, subd) ! get the owner process or subdomain of the element.
         if (DISTRIBUTED .and. (RANK .ne. subd)) cycle 
         select case(etype)

            case('mdlb')
               ! timer_a = MPI_Wtime()
               ! call project_p(mdle,flag_pref(iel),error_org,rate_p,Poly_flag)
               call project_p_linear(mdle,flag_pref(iel),error_org,rate_p,Poly_flag)
               ! timer_b = MPI_Wtime()
               call project_h(mdle,flag_pref(iel),error_org,rate_p,Poly_flag,Nref_grate(iel),Ref_indicator_flags((iel-1) * 10 + 1:iel*10))    
               ! timer_c = MPI_Wtime()
            case default
               write(*,*) "Element type not recognized"

         end select

         ! net_time_a = net_time_a + timer_b - timer_a
         ! net_time_b = net_time_b + timer_c - timer_b

      enddo
   !$OMP END PARALLEL DO 

   
   timer_save(istep,1) = net_time_a
   timer_save(istep,2) = net_time_b
   !MPI_REDUCTIONS
   count = nr_elem_ref

   call MPI_BARRIER (MPI_COMM_WORLD, ierr);

   call MPI_ALLREDUCE(MPI_IN_PLACE,Nref_grate,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   count = nr_elem_ref * 10
   call MPI_ALLREDUCE(MPI_IN_PLACE,Ref_indicator_flags,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

   ! write(*,*) "RANK = ",RANK," time for h and p projection = ",net_time_b,net_time_a
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   ! if(RANK .eq. ROOT) write(*,*) " time for h and p projection = ", end_time-start_time,net_time_a,net_time_b
   grate_mesh = maxval(Nref_grate)

   !---------------------------------------------------------------------
   if (RANK .eq. ROOT) write(*,*) "The guranteed rate = ",grate_mesh
   !---------------------------------------------------------------------------------------------
   ! copying back the old NODES Array 
   call Nodes_replace(NODES_cp)
   ! deallocate(NODES_cp)
!   allocate(NODES(MAXNODS))
   call move_alloc(NODES_cp, NODES)
   if(allocated(NODES_cp)) deallocate(NODES_cp)

   ! call dumpin_hp3d('befor_mesh_copy')
   NRELIS = NRELIS_old
   NRELES = NRELES_old
   NRNODS = NRNODS_old
   NRDOFSH = NRDOFSH_old
   NRDOFSE = NRDOFSE_old
   NRDOFSV = NRDOFSV_old
   NRDOFSQ = NRDOFSQ_old
   NPNODS = NPNODS_old
   MAXNODS = MAXNODS_old


   call update_ELEM_ORDER
   call enforce_min_rule
   call close_mesh
   call par_verify
   call update_gdof
   call update_Ddof


   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   !----------------------------------------------------------------------------------------------

   !! perfoming new refinements

   !--------------------------------------------------------------------------------
   select case(Irefine)

      case(IADAPTIVE)


         if (RANK .eq. ROOT)  write(*,*) 'Starting putting optimal p-refinement flags'
         if (RANK .eq. ROOT)  write(*,*) 'NRELES = ', NRELES  
      
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         !sets up p-refined flags for just p-refined elements
         ! allocate(write_pref(nr_elem_ref,2))
         ! write_pref = ZERO

         do iel = 1,nr_elem_ref

            mdle = mdle_ref(iel)
            etype = NODES(mdle)%type
            ! nord = NODES(mdle)%order
            select case(etype)

            case('mdlb')

               if(Nref_grate(iel) .gt. 0.25 * grate_mesh) then
                  if(Ref_indicator_flags((iel-1)*10 + 1) .eq. 2) then  !p-ref
                     
                     nord_new = Ref_indicator_flags((iel-1)*10 + 3)
                     write(*,*) "p-ref", mdle, nord_new
                     call nodmod(mdle,nord_new)
                     ! write_pref(iel,1) = mdle
                     ! write_pref(iel,2) = nord_new

                  endif
               endif
            
            case default
               write(*,*) "Element type not recognized"

            end select


         enddo


         ! open(4,file="p_ref_lst",status='replace')
         ! write(4,*) nr_elem_ref
         ! do iel = 1,nr_elem_ref
         !    write(4,*) write_pref(iel,:)
         ! enddo
         ! close(4)

         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         call enforce_min_rule
         call close_mesh
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call par_verify
         call update_gdof
         call update_Ddof


         call MPI_BARRIER (MPI_COMM_WORLD, ierr)

         ! putting h-refinement flags
         ! call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         if (RANK .eq. ROOT)  write(*,*) 'Starting to put  h-refined flags'
         ! allocate(write_kref(nr_elem_ref,2))
         ! write_kref = ZERO

            do iel = 1,nr_elem_ref

               mdle = mdle_ref(iel)
               etype = NODES(mdle)%type
               nord = NODES(mdle)%order
               select case(etype)

               case('mdlb')

                  if(Nref_grate(iel) .gt. 0.25 * grate_mesh) then
                     if(Ref_indicator_flags((iel-1)*10 + 1) .eq. 1) then   !h-ref

                        kref = Ref_indicator_flags((iel-1)*10 + 2)
                        if(RANK .eq. ROOT) write(*,*) "href = ",mdle,kref,NODES(mdle)%ref_kind
                        call refine(mdle,kref)
                        ! write_kref(iel,1) = mdle
                        ! write_kref(iel,2) = kref

                     endif
                  endif
               
               case default
                  write(*,*) "Element type not recognized"

               end select


            enddo



         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call close_mesh
         ! call enforce_min_rule
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call par_verify
         call update_gdof
         call update_Ddof
         call MPI_BARRIER (MPI_COMM_WORLD, ierr);end_time = MPI_Wtime()

         ! open(3,file="h_ref_lst",status='replace')
         ! write(3,*) nr_elem_ref
         ! do iel = 1,nr_elem_ref
         !    write(3,*) write_kref(iel,:)
         ! enddo
         ! close(3)

         if (RANK .eq. ROOT)  write(*,*) 'putting p-flags to the h-refined childs'

      
         ! !putting p-ref flags to the h-refined subchilds

         do iel = 1,nr_elem_ref

            mdle = mdle_ref(iel)
            etype = NODES(mdle)%type
            nord = NODES(mdle)%order
            first_son = NODES(mdle)%first_son

            select case(etype)

            case('mdlb')

               if(Nref_grate(iel) .gt. 0.25 * grate_mesh) then
                  if(Ref_indicator_flags((iel-1)*10 + 1) .eq. 1) then   !h-ref

                     kref_intent = Ref_indicator_flags((iel-1)*10+2)
                     kref_close  = NODES(mdle)%ref_kind

                     if(kref_intent .eq. kref_close) then

                        call ddecode(kref_intent,hx,hy,hz)
                        nr_sons_intent = 2**(hx+hy+hz)

                        do ic = 1,nr_sons_intent

                           nord_new = Ref_indicator_flags((iel-1)*10 + 2 + ic)
                           mdle_child = first_son + ic - 1
                           call nodmod(mdle_child,nord_new)

                        enddo
                     
                     else 

                        call ddecode(kref_intent,hx,hy,hz)
                        nr_sons_intent = 2**(hx+hy+hz)

                    
                        call ddecode(kref_close,hx,hy,hz)
                        nr_sons_close = 2**(hx+hy+hz)

                        allocate(pref_intent(nr_sons_intent))
                        allocate(pref_close(nr_sons_close))
                        pref_intent = ZERO
                        pref_close  = ZERO
                        ! write(*,*) nr_sons_intent, Ref_indicator_flags((iel-1)*10+3:(iel-1)*10+3+nr_sons_intent-1)
                        pref_intent(1:nr_sons_intent) = Ref_indicator_flags((iel-1)*10+3:(iel-1)*10+3+nr_sons_intent-1)
                     
                        call subson_one_irregularity_map(etype,kref_intent,kref_close,nr_sons_intent,pref_intent,nr_sons_close,pref_close)

                        do ic = 1,nr_sons_close

                           nord_new = pref_close(ic)
                           mdle_child = first_son + ic - 1
                           call nodmod(mdle_child,nord_new)

                        enddo

                        deallocate(pref_intent)
                        deallocate(pref_close)


                     endif

                  endif
               endif
            
            case default
               write(*,*) "Element type not recognized"

            end select



         enddo

         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call enforce_min_rule
         call close_mesh
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call par_verify
         call update_gdof
         call update_Ddof
         if (RANK .eq. ROOT)  write(*,*) 'end to redfinement'
         ! some printing of fnal flags for refined elements
         if (RANK .eq. ROOT) then
            do iel = 1,nr_elem_ref

               mdle = mdle_ref(iel)
               etype = NODES(mdle)%type
               nord = NODES(mdle)%order
               select case(etype)

               case('mdlb')

                  if(Nref_grate(iel) .gt. 0.25 * grate_mesh) then
                     ! if(Ref_indicator_flags(iel,1) .eq. 1) then   !h-ref

                        write(*,*) "The final ref flag = ",mdle,NODES(mdle)%ref_kind,NODES(mdle)%act
                     ! endif
                  endif
               
               case default
                  write(*,*) "Element type not recognized"

               end select


            enddo
            ! write(*,*)  "here = ",MAXNODS,NRNODS,NPNODS
         endif

         call MPI_BARRIER (MPI_COMM_WORLD, ierr)

      case default; Nstop = 1
   end select  
! endif
end subroutine HpAdapt





subroutine Hp_adapt_solve

   use common_prob_data
   use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
   use mpi_param     , only: RANK,ROOT
   use par_mesh      , only: DISTRIBUTED,HOST_MESH,distr_mesh
   use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_partition


   implicit none

   integer :: nsteps,count,src,ierr,i

   if (RANK .eq. ROOT) then
      nsteps=0
      do while (nsteps.le.0)
         write(*,*) 'Provide: number of adaptive hp-refinements'
         read(*,*) nsteps
      enddo
   else
      write(6,405) '[', RANK, '] : ','Waiting for broadcast from master...'
 405  format(A,I4,A,A)
   endif

   count = 1; src = ROOT
   call MPI_BCAST (nsteps,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)

      !..initial solve
   Print *, HOST_MESH, RANK
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      call par_mumps_sc('G')
   else
      call mumps_sc('G')
   endif

      !..initial solve
   Print *, HOST_MESH, RANK
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      call par_mumps_sc('G')
   else
      call mumps_sc('G')
   endif
   call exact_error
      ! ..do refinements and solve
   do i=0,nsteps
      !  ...display error and refine if necessary
            if (i.ne.nsteps) then
      !     ...adaptive refinement and solve
               call HpAdapt
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
               ! call refine_DPG
            endif
   enddo
   
end subroutine Hp_adapt_solve