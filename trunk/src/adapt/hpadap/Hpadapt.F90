!--------------------------------------------------------
!> @brief Routine performs a single step of Hp-refinement
!> @date June 2024
!--------------------------------------------------------
subroutine HpAdapt
!
   use control
   use data_structure3D
   use environment     , only: QUIET_MODE
   use assembly_sc     , only: NRDOF_CON, NRDOF_TOT
   use parametersDPG
   use par_mesh        , only: DISTRIBUTED, HOST_MESH, distr_mesh
   use mpi_param       , only: ROOT, RANK,NUM_PROCS
   use adaptivity      , only: CLOSURE_STRAT,MARKING_STRAT, &
                               ADAP_RATIO,ADAPT_STRAT
!
   implicit none
!
   integer, parameter :: max_step = 200
   integer, parameter :: IADAPTIVE = 2
!..Factor (when multiplied with the guaranteed rate)) that decides the
!  threshold for implementing anisotropic h or p or hp
   real(8), parameter :: Factor_max_rate = 0.25
!
   integer  :: Irefine = 2
   logical  :: ires = .true.
!
   integer :: Nflag(NR_PHYSA)
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
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
   real(8) :: resid_subd,resid_tot,elem_resid_max
   real(8) :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
   real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd
!
!..adaptive refinements
   integer :: nr_elem_ref
   integer :: ref_ndom(4)
   integer :: mdle_ref(NRELES)
   integer,allocatable  :: href_flags(:,:)
   real(8) :: res_sum
   integer,   allocatable :: ref_indicator_flags(:)
!..investment strategy variables
!..Currently works only for Hexahedral elements only
   integer,   allocatable :: mep_Nord_href(:,:,:)
   real(8),   allocatable :: error_reduction_rate(:,:)
   integer,   allocatable :: loc_max_rate(:)
   integer,   allocatable :: mep_count(:)
   real(8),   allocatable :: nref_grate(:)
   real(8) :: grate_mesh
   integer :: poly_flag
   integer,   allocatable :: flag_pref(:)
   integer :: nr_sons,first_son
   integer :: hx,hy,hz, nr_sons_intent,nr_sons_close,mdle_child
   integer, allocatable :: pref_intent(:), pref_close(:)
!
   integer :: i,ic,mdle,iel,kref,subd,count,ierr
   real(8) :: x(3), xnod(3,MAXbrickH)
   integer :: nord_new,is,nord
   real(8) :: error_org,rate_p
!
!..element type
   integer  :: etype
! for coarse mesh retrieval
   type(node), allocatable :: NODES_cp(:)
   integer :: NRELES_old,NRNODS_old,NRDOFSH_old,NRDOFSE_old,NRDOFSV_old,NRDOFSQ_old
   integer :: NPNODS_old,MAXNODS_old,NRELIS_old
   integer :: kref_close, kref_intent,href_count,counter
!
   real(8) :: start_time,end_time
!..printing flag
   integer :: iprint = 2
!..initialize
   elem_resid(1:NRELES)    = 0.d0
   elem_ref_flag(1:NRELES) = 0
   mdle_ref(1:NRELES)      = 0
   Nflag(1:NR_PHYSA) = (/0,0,1,1/) 
!
!..increase step if necessary
!..irefineold=0 means no refinement was made in the previous step
!..if first call or if a refinement was made, increase step
   if ((istep.eq.0).or.(irefineold.ne.0)) istep=istep+1
   irefineold=Irefine
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
   do iel = 1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
!
      if(ires) then
         call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
         resid_subd = resid_subd + elem_resid(iel) 
         if(elem_resid(iel) > elem_resid_max) elem_resid_max = elem_resid(iel)
      endif
!  
      if(NEXACT .ge. 1) then
         call element_error(mdle,Nflag,                     &
                        errorH,errorE,errorV,errorQ,       &
                        rnormH,rnormE,rnormV,rnormQ)
!
         error_subd = error_subd + errorQ
         rnorm_subd = rnorm_subd + rnormQ
!
      endif
!  
   enddo
!$OMP END PARALLEL DO
!
   resid_tot = 0.d0; error_tot = 0.d0; rnorm_tot = 0.d0
   if(DISTRIBUTED) then !MPI reduction op to collect the error and residual
      count = 1
      if(NEXACT.ge.1) then
         call MPI_ALLREDUCE(error_subd,error_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(rnorm_subd,rnorm_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif
      call MPI_ALLREDUCE(resid_subd,resid_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid_max,count,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr) !max over all elements
      count = NRELES
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   else
      if(NEXACT.ge.1) then
         error_tot = error_subd
         rnorm_tot = rnorm_subd
      endif
      resid_tot = resid_subd
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
   nelem_mesh(istep) =  NRELES
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
!
   open(22,file="error.dat",status='replace')
      do i  = 1,istep
         if(NEXACT .gt. 0) then
            write(22,*) i,nelem_mesh(i),nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i), &
                  error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
         else
            write(22,*) i,nelem_mesh(i),nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i)
         endif
      enddo
   close(22)
!
   70 continue
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
!..mark elements for adaptive refinement
   if (Irefine .eq. IADAPTIVE) then
      nr_elem_ref = 0
!  ...strategy 1 (greedy strategy)
      if(MARKING_STRAT .eq. 0) then
         do iel = 1,NRELES
               mdle = ELEM_ORDER(iel)
               call find_domain(mdle,i)
               if(elem_resid(iel) > ADAP_RATIO * elem_resid_max) then
                  nr_elem_ref = nr_elem_ref + 1
                  mdle_ref(nr_elem_ref) = mdle
                  if (RANK.eq.ROOT) write(*,1020) 'ndom=',i,', mdle=',mdle,', r =',elem_resid(iel)
               endif
               1020    format(A,I2,A,I9,A,es12.5)
         enddo
      elseif(MARKING_STRAT .eq. 1) then    ! strategy 2: Doerfler's Strategy
!    
         mdle_ref(1:NRELES) = ELEM_ORDER(1:NRELES)
         call qsort_duplet(mdle_ref,elem_resid,NRELES,1,NRELES)
!..verify the array is sorted
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
!    
         res_sum = 0.d0
         ref_ndom(1:4)=0
!    
         do iel = 1,NRELES
               mdle = mdle_ref(iel)
               call find_domain(mdle,i)
               ref_ndom(i) = ref_ndom(i) + 1
!    
               nr_elem_ref = nr_elem_ref + 1
               res_sum = res_sum + elem_resid(iel)
               if(res_sum > ADAP_RATIO*resid_tot) exit
!    
         enddo
      endif
   endif
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!..start copying the current nodes structure into a copy nodes structure
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
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!..flag_pref contains the flags if element has been hp refined
   allocate(flag_pref(nr_elem_ref))
   flag_pref = ZERO
!..hp refining the marked elements
!..setting up the ANSIO_ADAP flag from adaptation module to allow isotropic adaptation
   CLOSURE_STRAT = 0
   call Finegrid_padap(nr_elem_ref,mdle_ref,flag_pref)
!
!..solving on the hp refined mesh to fine mesh solution
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      call par_mumps_sc('G')
   else
      call mumps_sc('G')
   endif
!..exact error on the Hp refined reference mesh
   if(NEXACT .eq. 1) call exact_error
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..allocating memory for arrays required in Hp adaptation
   allocate(mep_Nord_href(100,8,nr_elem_ref))
   mep_Nord_href = ZERO
!  
   allocate(error_reduction_rate(100,nr_elem_ref))
   error_reduction_rate = ZERO
!
   allocate(loc_max_rate(nr_elem_ref))
   loc_max_rate = ZERO
!
   allocate(mep_count(nr_elem_ref))
   mep_count = ZERO
!
!..array to store href flags for elements which are selected for h-refinement for use in refine_opt
   allocate(href_flags(nr_elem_ref,2))
   href_flags = ZERO
!
   allocate(nref_grate(nr_elem_ref))
   nref_grate = ZERO
!
!..Note regarding ref_indicator_flags array: first column will contain indicator flag for
!  kind of refinement: 1 for h-ref and 2 for p-ref, 2nd column will have the refinement flag: eg: 000 means no
!  h-ref and only p-ref and non-zero href indicator means href, from 3rd column onwards we have the polynomial
!  distribution for each subson generated during h refinement, Note: the memory allocation for p-dist of subson
!  is done for isotropic brick element as it will produce maximum subsons and thus, memory alocation will work for
!  every other element like prism, pyramid and tets.
!
   allocate(ref_indicator_flags(nr_elem_ref * 10))
   ref_indicator_flags = ZERO
!
!$OMP PARALLEL DO                                          &
!$OMP PRIVATE(mdle,etype,subd,error_org,rate_p,poly_flag)  &
!$OMP SCHEDULE(DYNAMIC)
   do iel = 1,nr_elem_ref
      mdle = mdle_ref(iel)
      etype = NODES(mdle)%ntype
      call get_subd(mdle, subd) ! get the owner process or subdomain of the element.
      if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
      select case(etype)
!
      case(MDLB)
         call project_p_linear(mdle,flag_pref(iel),error_org,rate_p,poly_flag,istep)
!
         if(ADAPT_STRAT .eq. 1) then
            call project_h(mdle,flag_pref(iel),error_org,rate_p,poly_flag,istep, &
                           nref_grate(iel),ref_indicator_flags((iel-1) * 10 + 1:iel*10), &
                           mep_Nord_href(1:100,1:8,iel),error_reduction_rate(1:100,iel), &
                           mep_count(iel),loc_max_rate(iel))
         endif
!
      case default
         write(*,*) "Element type not recognized"
      end select
   enddo
!$OMP END PARALLEL DO 
!
!..MPI_REDUCTIONS
   count = nr_elem_ref
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,nref_grate,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   grate_mesh = maxval(nref_grate)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr);
!
!..loop for deciding the polynomial distribution in case of h-refinements using investment strat
!$OMP PARALLEL DO                                          &
!$OMP PRIVATE(mdle,etype,subd,error_org,rate_p,poly_flag) &
!$OMP SCHEDULE(DYNAMIC)
   do iel = 1,nr_elem_ref  
      mdle = mdle_ref(iel)
      etype = NODES(mdle)%ntype        
      call get_subd(mdle, subd) ! get the owner process or subdomain of the element.
      if (DISTRIBUTED .and. (RANK .ne. subd)) cycle 
      select case(etype)
!
      case(MDLB)
         if(ref_indicator_flags((iel-1)*10 + 1) .eq. 1) then
            ref_indicator_flags((iel-1)*10+3:iel*10) = ZERO               
            call dof_investment(Factor_max_rate, grate_mesh,mep_Nord_href(1:100,1:8,iel), &
                                 error_reduction_rate(1:100,iel),loc_max_rate(iel), &
                                 mep_count(iel),ref_indicator_flags((iel-1)*10+1:iel*10))
         endif
      case default
         write(*,*) "Element type not recognized"
!
      end select
   enddo
!$OMP END PARALLEL DO
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr);
   count = nr_elem_ref * 10
   call MPI_ALLREDUCE(MPI_IN_PLACE,ref_indicator_flags,count,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!
!---------------------------------------------------------------------
   if (RANK .eq. ROOT) write(*,*) "The guranteed rate = ",grate_mesh
!---------------------------------------------------------------------
!..copying back the old NODES Array 
   call Nodes_dealloc(NODES_cp)
   call move_alloc(NODES_cp, NODES)
   if(allocated(NODES_cp)) deallocate(NODES_cp)
!
   NRELIS = NRELIS_old
   NRELES = NRELES_old
   NRNODS = NRNODS_old
   NRDOFSH = NRDOFSH_old
   NRDOFSE = NRDOFSE_old
   NRDOFSV = NRDOFSV_old
   NRDOFSQ = NRDOFSQ_old
   NPNODS = NPNODS_old
   MAXNODS = MAXNODS_old
!
   call update_ELEM_ORDER
   call enforce_min_rule
   call close_mesh
   call par_verify
   call update_gdof
   call update_Ddof
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!..setting up the ANSIO_ADAP flag from adaptation module to allow aniso adaptation
   CLOSURE_STRAT = 1
!----------------------------------------------------------------------------------------------
!..perfoming computed Hp refinements
!--------------------------------------------------------------------------------
   select case(Irefine)
      case(IADAPTIVE)
!
         if (RANK .eq. ROOT)  write(*,*) 'Starting putting optimal p-refinement flags'
         if (RANK .eq. ROOT)  write(*,*) 'NRELES = ', NRELES  
!        
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!     ...sets up p-refined flags for just p-refined elements
         href_count = 0
         do iel = 1,nr_elem_ref
!
            mdle = mdle_ref(iel)
            etype = NODES(mdle)%ntype
            select case(etype)
            case(MDLB)
               if(nref_grate(iel) .gt. Factor_max_rate * grate_mesh) then
                  if(ref_indicator_flags((iel-1)*10 + 1) .eq. 2) then
!                        
                        nord_new = ref_indicator_flags((iel-1)*10 + 3)
                        call nodmod(mdle,nord_new)
                        href_count = href_count + 1
                  endif
               endif
            case default
               write(*,*) "Element type not recognized"
            end select
         enddo
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         call enforce_min_rule
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call par_verify
         call update_gdof
         call update_Ddof
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!     ...computing the href_flags list
         counter = 1
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
            etype = NODES(mdle)%ntype
            nord = NODES(mdle)%order
            select case(etype)
            case(MDLB)
               if(nref_grate(iel) .gt. Factor_max_rate * grate_mesh) then
                  if(ref_indicator_flags((iel-1)*10 + 1) .eq. 1) then
!
                     kref = ref_indicator_flags((iel-1)*10 + 2)
                     href_flags(counter,1) = mdle
                     href_flags(counter,2) = kref
                     counter = counter + 1
                  endif
               endif
            case default
               write(*,*) "Element type not recognized"
            end select
         enddo
         counter = counter - 1
!     ...putting h-refinement flags
         if (RANK .eq. ROOT)  write(*,*) 'Starting to put  h-refined flags'
         href_count = 0
!
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
            etype = NODES(mdle)%ntype
            nord = NODES(mdle)%order
            select case(etype)
!
            case(MDLB)
               if(nref_grate(iel) .gt. Factor_max_rate * grate_mesh) then
                  if(ref_indicator_flags((iel-1)*10 + 1) .eq. 1) then   !h-ref
!
                     kref = ref_indicator_flags((iel-1)*10 + 2)
                     call refine_opt(mdle,kref,href_flags(1:counter,:),counter)
                     href_count = href_count + 1
!
                  endif
               endif
            case default
               write(*,*) "Element type not recognized"
               stop
            end select
         enddo
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call close_mesh
         call enforce_min_rule
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call par_verify
         call update_gdof
         call update_Ddof
         call MPI_BARRIER (MPI_COMM_WORLD, ierr);end_time = MPI_Wtime()
!
         if (RANK .eq. ROOT)  write(*,*) 'putting p-flags to the h-refined childs'
!     ...putting p-ref flags to the h-refined subchilds
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
            etype = NODES(mdle)%ntype
            nord = NODES(mdle)%order
            first_son = NODES(mdle)%first_son
            select case(etype)
            case(MDLB)
               if(nref_grate(iel) .gt. Factor_max_rate * grate_mesh) then
                  if(ref_indicator_flags((iel-1)*10 + 1) .eq. 1) then   !h-ref
!
                     kref_intent = ref_indicator_flags((iel-1)*10+2)
                     kref_close  = NODES(mdle)%ref_kind
                     if(kref_intent .eq. kref_close) then
                        call ddecode(kref_intent,hx,hy,hz)
                        nr_sons_intent = 2**(hx+hy+hz)
!
                        do ic = 1,nr_sons_intent
                           nord_new = ref_indicator_flags((iel-1)*10 + 2 + ic)
                           mdle_child = first_son + ic - 1
                           call nodmod(mdle_child,nord_new)
                        enddo
!
                     else 
                        call ddecode(kref_intent,hx,hy,hz)
                        nr_sons_intent = 2**(hx+hy+hz)
                        call ddecode(kref_close,hx,hy,hz)
                        nr_sons_close = 2**(hx+hy+hz)
                        allocate(pref_intent(nr_sons_intent))
                        allocate(pref_close(nr_sons_close))
                        pref_intent = ZERO
                        pref_close  = ZERO
                        pref_intent(1:nr_sons_intent) = ref_indicator_flags((iel-1)*10+3:(iel-1)*10+3+nr_sons_intent-1)
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
               stop
            end select
         enddo
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call close_mesh
         call enforce_min_rule
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call par_verify
         call update_gdof
         call update_Ddof
   end select  
!..deallocating the allocatable arrays
   deallocate(href_flags)
   deallocate(ref_indicator_flags)
   deallocate(mep_Nord_href)
   deallocate(error_reduction_rate)
   deallocate(loc_max_rate)
   deallocate(mep_count)
   deallocate(nref_grate)
   deallocate(flag_pref)
!
end subroutine HpAdapt
!----------------------------------------------------------------------
!> @brief Driver routine to run Hp adaptation.
!> @date June 2024
!----------------------------------------------------------------------
subroutine Hp_adapt_solve
!
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
         write(*,*) 'Provide: number of adaptive hp-refinements'
         read(*,*) nsteps
      enddo
   else
      write(6,405) '[', RANK, '] : ','Waiting for broadcast from master...'
   405  format(A,I4,A,A)
   endif
!
   count = 1; src = ROOT
   call MPI_BCAST (nsteps,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)

!..initial solve
   Print *, HOST_MESH, RANK
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      call par_mumps_sc('G')
   else
      call mumps_sc('G')
   endif
!..do refinements and solve
   do i=0,nsteps
!..display error and refine if necessary
      if (i.ne.nsteps) then
!..adaptive refinement and solve
         call HpAdapt
         if (DISTRIBUTED .and. (.not. HOST_MESH)) then
               call par_mumps_sc('G')
         else
               call mumps_sc('G')
         endif
      endif
   enddo
!
end subroutine Hp_adapt_solve
