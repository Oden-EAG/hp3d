!-----------------------------------------------------------------------------
!> @brief Enforces one-irregular mesh
!> @date Feb 2023
subroutine close_mesh()
   use error
   use refinements
   use data_structure3D
   use par_mesh   , only: DISTRIBUTED
   use environment, only: QUIET_MODE
   use mpi_wrapper
   use bitvisit
!
   implicit none
!
   integer                :: ntype
   integer, allocatable   :: list(:,:)
   integer, dimension(27) :: nodesl,norientl
   integer, dimension(6)  :: kreff
   integer, dimension(12) :: krefe
   integer :: i, j, ic, mdle, nod, kref
   integer :: nreles_aux
   logical :: nflag
   real(8) :: start_time,end_time
   integer :: ierr
   integer :: nv,ne,nf,nve
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!-----------------------------------------------------------------------------
!
   if (DISTRIBUTED) then
      call close_mesh_par
      goto 99
   endif
!
   do
!---------------------------------------------------------
!     Step 0 : activate nodes
!---------------------------------------------------------
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
      call refresh
      call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
      if (.not.QUIET_MODE .and. RANK.eq.ROOT) then
         write(*,2020) end_time-start_time
    2020 format(' refresh    : ',f12.5,'  seconds')
      endif
!
!---------------------------------------------------------
!     Step 1 : check constrained nodes
!---------------------------------------------------------
#if HP3D_DEBUG
      if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
         write(*,*) 'close_mesh: CHECK CONSTRAINING NODES'
      endif
#endif
!
      call bitvisit_init(NRNODS)
!
!  ...allocate list
      allocate(list(2,NRELES))
      list(1:2,1:NRELES) = 0; ic = 0
!
!$OMP PARALLEL
!$OMP DO PRIVATE(mdle,nodesl,norientl)
      do i=1,NRELES
         mdle = ELEM_ORDER(i)
         call get_connect_info(mdle, nodesl,norientl) ! setting internal arrays
         call flag_constr_parents(mdle,nodesl)
      enddo
!$OMP END DO
!
#if HP3D_DEBUG
   if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
      !$OMP SINGLE
      write(*,*) 'close_mesh: RESOLVE THE NODES MORE THAN ONE IRREGULARITY'
      !$OMP END SINGLE
   endif
#endif
!
!---------------------------------------------------------
!     Step 2 : pick the middle nodes to be refined
!---------------------------------------------------------
!
!$OMP DO                                              &
!$OMP PRIVATE(mdle,nodesl,norientl,nod,ntype,nflag,   &
!$OMP         krefe,kreff,kref,j,nv,ne,nf,nve)        &
!$OMP REDUCTION (+:ic)
      do i=1,NRELES
         mdle = ELEM_ORDER(i)
         call elem_nodes(mdle, nodesl,norientl)
!
         ntype = NODES(mdle)%ntype
         nflag = .false.
!
!        check edges
!        ~~~~~~~~~~~~
         krefe = 0
!
         nv = nvert(ntype)
         ne = nedge(ntype)
         nf = nface(ntype)
         nve = nv + ne
!
         do j = 1,ne
            nod = nodesl(nv+j)
            if (visited(nod)) then
               krefe(j)=1
               nflag=.true.
            endif
         enddo
!
!        check faces
!        ~~~~~~~~~~~~
         kreff = 0
         do j = 1,nf
            nod = nodesl(nve+j)
            if (visited(nod)) then
               call get_isoref(nod, kreff(j))
               nflag=.true.
            endif
         enddo
!
!        find out more than one-irregular
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (nflag) then
!
#if HP3D_DEBUG
            if ((RANK.eq.ROOT) .and. (iprint.eq.1)) then
               !$OMP CRITICAL
               write(*,*) 'close_mesh: mdle = ', mdle
               write(*,7003) krefe(1:12)
   7003        format('krefe = ',12i2)
               write(*,7004) kreff(1:6)
   7004        format('kreff = ',6i3)
               call pause
               !$OMP END CRITICAL
            endif
#endif
!           increment counter
            ic=ic+1
!
!           find out proper refinement flag
            if (is_iso_only()) then
               call get_isoref(mdle, kref)
            else
!           ...-------------------------------------------------------------------
!           ...Option 1: do minimum refinement that is necessary
!               call find_element_closing_ref(ntype,kreff,krefe, kref)
!           ...-------------------------------------------------------------------
!           ...Option 2: always ask for isotropic refinement
               call get_isoref(mdle, kref)
!           ...-------------------------------------------------------------------
!           ...Option 3: always ask for radial (xy) refinement (FIBER LASER)
!               select case (NODES(mdle)%ntype)
!                  case(MDLB); kref = 110
!                  case(MDLP); kref = 10
!               end select
!           ...-------------------------------------------------------------------
            endif
            list(1,i) = mdle
            list(2,i) = kref
         endif
      enddo
!$OMP END DO
!$OMP END PARALLEL
!
      call bitvisit_finalize
!
#if HP3D_DEBUG
      if ((.not.QUIET_MODE.or.iprint.eq.2) .and. (RANK.eq.ROOT)) then
         write(*,7002) ' close_mesh : number of elements to refine = ', ic
   7002  format(A,i6)
      endif
#endif
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
#if HP3D_DEBUG
            if ((RANK.eq.ROOT) .and. (iprint.eq.1)) then
               write(*,7001) i, mdle, S_Type(NODES(mdle)%ntype), kref
    7001       format('close_mesh: i= ',i6,' mdle= ', i6,' ', a4, ' ref_kind = ',i5)
            endif
#endif
            call refine(mdle,kref)
         endif
      enddo
      deallocate(list)
!..end outer loop
   enddo
!
   if (allocated(list)) deallocate(list)
!
   99 continue
!
!..ensure that DOFs are allocated correctly within subdomains
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   call distr_refresh
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (.not.QUIET_MODE .and. RANK.eq.ROOT) then
      write(*,2030) end_time-start_time
 2030 format(' distr_refr : ',f12.5,'  seconds')
   endif
!
#if HP3D_DEBUG
   call par_verify
#endif
!
!
end subroutine close_mesh
!
!
!
!-----------------------------------------------------------------------------
!> @brief enforce one-irregular mesh
!> @note MPI-distributed version of close_mesh
!> @date Mar 2023
subroutine close_mesh_par()
   use error
   use refinements
   use data_structure3D
   use mpi_wrapper
   use bitvisit
   use environment, only: QUIET_MODE
!
   implicit none
!
   integer                :: ntype
   integer, allocatable   :: list(:,:), list_glob(:,:)
   integer, dimension(27) :: nodesl,norientl
   integer, dimension(6)  :: kreff
   integer, dimension(12) :: krefe
!
   integer :: i, j, ic, mdle, nod, kref, ierr
   integer :: nv,ne,nf,nve
   integer :: ic_glob
   logical :: nflag
   real(8) :: start_time,end_time
   integer :: displs(NUM_PROCS), ic_procs(NUM_PROCS)
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
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
      if (.not.QUIET_MODE .and. RANK.eq.ROOT) write(*,2020) end_time-start_time
 2020 format(' refresh    : ',f12.5,'  seconds')
!
!---------------------------------------------------------
!     Step 1 : check constrained nodes
!---------------------------------------------------------
#if HP3D_DEBUG
      if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
         write(*,*) 'close_mesh_par: CHECK CONSTRAINING NODES'
      endif
#endif
!
!  ...bit flag is like NODES(nod)%visit but is stored as bit collection
      call bitvisit_init(NRNODS)
!
!  ...allocate list
      allocate(list(2,NRELES_SUBD))
      list(1:2,1:NRELES_SUBD) = 0; ic = 0
!
!$OMP PARALLEL DO PRIVATE(mdle,nodesl,norientl)
      do i=1,NRELES_SUBD
         mdle = ELEM_SUBD(i)
         call get_connect_info(mdle, nodesl,norientl) ! setting internal arrays
         call flag_constr_parents(mdle,nodesl)
      enddo
!$OMP END PARALLEL DO
!
!  ...reduce flagged nodes
!  ...Note: could instead loop over NRELES_GHOST but this is more simple
!           for now.
      call reduce_visit
!
#if HP3D_DEBUG
   if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
      write(*,*) 'close_mesh_par: RESOLVE THE NODES MORE THAN ONE IRREGULARITY'
   endif
#endif
!
!---------------------------------------------------------
!     Step 2 : pick the middle nodes to be refined
!---------------------------------------------------------
!
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(mdle,nodesl,norientl,nod,ntype,nflag,   &
!$OMP         krefe,kreff,kref,j,nv,ne,nf,nve)
      do i=1,NRELES_SUBD
         mdle=ELEM_SUBD(i)
         call elem_nodes(mdle, nodesl,norientl)
!
         ntype = NODES(mdle)%ntype
         nflag = .false.
!
!        check edges
!        ~~~~~~~~~~~~
         krefe = 0
!
         nv = nvert(ntype)
         ne = nedge(ntype)
         nf = nface(ntype)
         nve = nv + ne
!
         do j = 1,ne
            nod = nodesl(nv+j)
            if (visited(nod)) then
               krefe(j)=1
               nflag=.true.
            endif
         enddo
!
!        check faces
!        ~~~~~~~~~~~~
         kreff = 0
         do j = 1,nf
            nod = nodesl(nve+j)
            if (visited(nod)) then
               call get_isoref(nod, kreff(j))
               nflag=.true.
            endif
         enddo
!
!        find out more than one-irregular
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (nflag) then
!
#if HP3D_DEBUG
            if ((RANK.eq.ROOT) .and. (iprint.eq.1)) then
               !$OMP CRITICAL
               write(*,*) 'close_mesh_par: mdle = ', mdle
               write(*,7003) krefe(1:12)
   7003        format('krefe = ',12i2)
               write(*,7004) kreff(1:6)
   7004        format('kreff = ',6i3)
               call pause
               !$OMP END CRITICAL
            endif
#endif
!
!           find out proper refinement flag
            if (is_iso_only()) then
               call get_isoref(mdle, kref)
            else
!           ...-------------------------------------------------------------------
!           ...Option 1: do minimum refinement that is necessary
!               call find_element_closing_ref(ntype,kreff,krefe, kref)
!           ...-------------------------------------------------------------------
!           ...Option 2: always ask for isotropic refinement
               call get_isoref(mdle, kref)
!           ...-------------------------------------------------------------------
!           ...Option 3: always ask for radial (xy) refinement (FIBER LASER)
!               select case (NODES(mdle)%ntype)
!                  case(MDLB); kref = 110
!                  case(MDLP); kref = 10
!               end select
!           ...-------------------------------------------------------------------
            endif
!
!$OMP CRITICAL
!           increment counter
            ic=ic+1
!
            list(1,ic) = mdle
            list(2,ic) = kref
!$OMP END CRITICAL
         endif
      enddo
!$OMP END PARALLEL DO
!
      call bitvisit_finalize
!
      call MPI_Allgather(ic,1,MPI_INTEGER,ic_procs,1,MPI_INTEGER,MPI_COMM_WORLD, ierr)
!
      displs = 0
      ic_glob = 0
      do i=1,NUM_PROCS
         displs(i) = 2*ic_glob
         ic_glob = ic_glob + ic_procs(i)
      enddo
!
#if HP3D_DEBUG
      if (RANK.eq.ROOT .and. iprint.eq.2) then
         write(*,*) 'close_mesh_par: number of elements to refine ', ic_glob
      endif
#endif
!
!     loop exit condition
      if (ic_glob.eq.0) exit
!
      allocate(list_glob(2,ic_glob))
      call MPI_AllgatherV(list(1:2,1:ic),2*ic,MPI_INTEGER,list_glob, 2*ic_procs,displs,MPI_INTEGER,MPI_COMM_WORLD, ierr)
      deallocate(list)
!
!---------------------------------------------------------
!     Step : refine from the list
!---------------------------------------------------------
      do i=1,ic_glob
         mdle = list_glob(1,i)
         kref = list_glob(2,i)
         if (is_leaf(mdle)) then
#if HP3D_DEBUG
            if ((RANK.eq.ROOT) .and. (iprint.eq.1)) then
               write(*,7001) i, mdle, S_Type(NODES(mdle)%ntype), kref
    7001       format('close_mesh_par: i= ',i6,' mdle= ', i6,' ', a4, ' ref_kind = ',i5)
            endif
#endif
            call refine(mdle,kref)
         endif
      enddo
      deallocate(list_glob)
!..end outer loop
   enddo
!
   if (allocated(list)) deallocate(list)
!
end subroutine close_mesh_par
