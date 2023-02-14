

! routines for local checks: will not be pushed in


subroutine close_mesh_check()

    use error
    use refinements
    use data_structure3D
    use mpi_param  , only: ROOT,RANK
    use MPI        , only: MPI_COMM_WORLD
 !
    implicit none
 !
    character(4)           :: type
    integer, allocatable   :: list(:,:)
    integer, dimension(27) :: nodesl,norientl
    integer, dimension(6)  :: kreff
    integer, dimension(12) :: krefe
    integer :: i, j, ic, mdle, nod, kref
    integer :: nreles_aux
    logical :: nflag
    real(8) :: MPI_Wtime,start_time,end_time
    integer :: ierr
 !
!  #if DEBUG_MODE
!     integer :: nre, nrf
!     integer :: iprint = 0
!  #endif
 !
 !-----------------------------------------------------------------------------
 !
    do
 !---------------------------------------------------------
 !     Step 0 : activate nodes
 !---------------------------------------------------------
       call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
       call refresh
       call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
       if (RANK .eq. ROOT) write(*,2020) end_time-start_time
  2020 format(' refresh    : ',f12.5,'  seconds')
 !
 !---------------------------------------------------------
 !     Step 1 : check constrained nodes
 !---------------------------------------------------------
!  #if DEBUG_MODE
      !  if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
      !     write(*,*) 'close_mesh: CHECK CONSTRAINING NODES'
      !  endif
!  #endif
 !
 !  ...allocate list
       allocate(list(2,NRELES))
       list(1:2,1:NRELES) = 0; ic = 0
 !
 !$OMP PARALLEL
 !$OMP DO
       do i=1,NRNODS
         NODES(i)%visit = 0
       enddo
 !$OMP END DO
 !$OMP DO PRIVATE(mdle,nodesl,norientl)
       do i=1,NRELES
          mdle = ELEM_ORDER(i)
          call get_connect_info(mdle, nodesl,norientl) ! setting internal arrays
          call flag_constr_parents(mdle)
       enddo
 !$OMP END DO
 !
!  #if DEBUG_MODE
   !  if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
   !     !$OMP SINGLE
   !     write(*,*) 'close_mesh: RESOLVE THE NODES MORE THAN ONE IRREGULARITY'
   !     !$OMP END SINGLE
   !  endif
!  #endif
 !
 !---------------------------------------------------------
 !     Step 2 : pick the middle nodes to be refined
 !---------------------------------------------------------
 !
 !$OMP DO                                     &
 !$OMP PRIVATE(mdle,nodesl,norientl,nod,type, &
 !$OMP         nflag,krefe,kreff,kref,j)      &
 !$OMP REDUCTION (+:ic)
       do i=1,NRELES
          mdle=ELEM_ORDER(i)
          call elem_nodes(mdle, nodesl,norientl)
 !
          type  = NODES(mdle)%type
          nflag = .FALSE.
 !
 !        check edges
 !        ~~~~~~~~~~~~
          krefe = 0
          do j=1, nedge(type)
             nod = nodesl(nvert(type)+j)
             if (NODES(nod)%visit.eq.1) then
                krefe(j)=1
                nflag=.TRUE.
             endif
          enddo
 !
 !        check faces
 !        ~~~~~~~~~~~~
          kreff = 0
          do j=1, nface(type)
             nod = nodesl(nvert(type)+nedge(type)+j)
             if (NODES(nod)%visit.eq.1) then
                call get_isoref(nod, kreff(j))
                nflag=.TRUE.
             endif
          enddo
 !
 !        find out more than one-irregular
 !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (nflag) then
 !
!  #if DEBUG_MODE
!              if ((RANK.eq.ROOT) .and. (iprint.eq.1)) then
!                 !$OMP CRITICAL
!                 nre = nedge(NODES(mdle)%type)
!                 nrf = nface(NODES(mdle)%type)
!                 write(*,*) 'close_mesh: mdle = ', mdle
!                 write(*,7003) krefe(1:12)
!     7003        format('krefe = ',12i2)
!                 write(*,7004) kreff(1:6)
!     7004        format('kreff = ',6i3)
!                 call pause
!                 !$OMP END CRITICAL
!              endif
!  #endif
 !           increment counter
             ic=ic+1
 !
 !           find out proper refinement flag
             if (is_iso_only()) then
                call get_isoref(mdle, kref)
             else
 !           ...-------------------------------------------------------------------
 !           ...Option 1: do minimum refinement that is necessary
                !call find_element_closing_ref(type,kreff,krefe, kref)
 !           ...-------------------------------------------------------------------
 !           ...Option 2: always ask for isotropic refinement
                call get_isoref(mdle, kref)
 !           ...-------------------------------------------------------------------
 !           ...Option 3: always ask for radial (xy) refinement (FIBER LASER)
                !select case (NODES(mdle)%type)
                !   case('mdlb'); kref = 110
                !   case('mdlp'); kref = 10
                !end select
 !           ...-------------------------------------------------------------------
             endif
             list(1,i) = mdle
             list(2,i) = kref
          endif
       enddo
 !$OMP END DO
 !$OMP END PARALLEL
 !
!  #if DEBUG_MODE
!        iprint = 2
!        if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
!           write(*,*) 'close_mesh: number of elements to refine ', ic
!        endif
!        iprint = 0
!  #endif
 !
 !     loop exit condition
       if (ic.eq.0) exit
 !
 !---------------------------------------------------------
 !     Step : refine from the list
 !---------------------------------------------------------
       nreles_aux = NRELES
       do i=1,nreles_aux
          if (list(1,i) .eq. 0) cycle
          mdle = list(1,i)
          kref = list(2,i)
          if (is_leaf(mdle)) then
!  #if DEBUG_MODE
!              if ((RANK.eq.ROOT) .and. (iprint.eq.1)) then
!                 write(*,7001) i, mdle, NODES(mdle)%type, kref
!      7001       format('close_mesh: i= ',i6,' mdle= ', i6,' ', a4, ' ref_kind = ',i5)
!              endif
!  #endif
             call refine(mdle,kref)
          endif
       enddo
       deallocate(list)
 !..end outer loop
    enddo
 !
    if (allocated(list)) deallocate(list)

end subroutine close_mesh_check


subroutine restart_adaptation()

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

   
   if (RANK .eq. ROOT) call dumpin_hp3d("10adp_mesh_copy")
   call update_ELEM_ORDER
   call enforce_min_rule
   call par_verify
   call close_mesh
   call update_gdof
   call update_Ddof


end subroutine restart_adaptation


subroutine restart_adaptation_read_ref()

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
   character(len=12) :: Dump_file
   integer :: ndmp
   integer :: nr_elem_ref,iel,val1,val2,kref,mdle,pref
   integer, allocatable :: href_list(:,:),pref_list(:,:)
   Dump_file = "adap_ref_lst"
   
   if (RANK .eq. ROOT) call dumpin_hp3d("10adp_mesh_copy")
   call update_ELEM_ORDER
   call enforce_min_rule
   call par_verify
   call close_mesh
   call update_gdof
   call update_Ddof

   !read list
   ndmp = 3
   open(unit = ndmp,file = Dump_file)
   read(ndmp,*) nr_elem_ref
   allocate(href_list(nr_elem_ref,2))

   do iel = 1,nr_elem_ref
      read(ndmp,*) val1,val2
      href_list(iel,1) = val1
      href_list(iel,2) = val2
      write(*,*) href_list(iel,:)
   enddo
   
   close(3)

   ndmp = 4
   Dump_file = "p_ref_lst"
   open(unit = ndmp,file = Dump_file)
   read(ndmp,*) nr_elem_ref
   allocate(pref_list(nr_elem_ref,2))

   do iel = 1,nr_elem_ref
      read(ndmp,*) val1,val2
      pref_list(iel,1) = val1
      pref_list(iel,2) = val2
      write(*,*) href_list(iel,:)
   enddo
   
   close(4)

   do iel = 1,nr_elem_ref
      if(pref_list(iel,1) .gt. 0) then
         mdle = pref_list(iel,1)
         pref = pref_list(iel,2)
         call nodmod(mdle,pref)
      endif
   enddo

   call enforce_min_rule
   call close_mesh
   ! call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   call par_verify
   call update_gdof
   call update_Ddof

   do iel = 1,3
      if(href_list(iel,1) .gt. 0) then
         kref = href_list(iel,2)
         mdle = href_list(iel,1)
         call refine(mdle,kref)
      endif
   enddo
   ! call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   call close_mesh
   call enforce_min_rule
   ! call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   call par_verify
   call update_gdof
   call update_Ddof

end subroutine restart_adaptation_read_ref