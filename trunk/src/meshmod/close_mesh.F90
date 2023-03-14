subroutine close()
   implicit none
   integer :: iprint
!..ensure backward compatibility
   iprint=0
   if (iprint /= 0)  write(*,*) 'Mesh closing begin'
   call close_mesh
   if (iprint /= 0)  write(*,*) 'Mesh closing end'
end subroutine close
!
!-----------------------------------------------------------------------------
!> @brief enforce one-irregular mesh
!> @date Feb 2023
subroutine close_mesh()
   use error
   use refinements
   use data_structure3D
   use mpi_param  , only: ROOT,RANK
   use MPI        , only: MPI_COMM_WORLD
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
   real(8) :: MPI_Wtime,start_time,end_time
   integer :: ierr
!
#if DEBUG_MODE
   integer :: nre, nrf
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
      if (RANK .eq. ROOT) write(*,2020) end_time-start_time
 2020 format(' refresh    : ',f12.5,'  seconds')
!
!---------------------------------------------------------
!     Step 1 : check constrained nodes
!---------------------------------------------------------
#if DEBUG_MODE
      if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
         write(*,*) 'close_mesh: CHECK CONSTRAINING NODES'
      endif
#endif
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
         call flag_constr_parents(mdle,nodesl)
      enddo
!$OMP END DO
!
#if DEBUG_MODE
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
!$OMP DO                                      &
!$OMP PRIVATE(mdle,nodesl,norientl,nod,ntype, &
!$OMP         nflag,krefe,kreff,kref,j)       &
!$OMP REDUCTION (+:ic)
      do i=1,NRELES
         mdle=ELEM_ORDER(i)
         call elem_nodes(mdle, nodesl,norientl)
!
         ntype = NODES(mdle)%ntype
         nflag = .FALSE.
!
!        check edges
!        ~~~~~~~~~~~~
         krefe = 0
         do j=1,nedge(ntype)
            nod = nodesl(nvert(ntype)+j)
            if (NODES(nod)%visit.eq.1) then
               krefe(j)=1
               nflag=.TRUE.
            endif
         enddo
!
!        check faces
!        ~~~~~~~~~~~~
         kreff = 0
         do j=1,nface(ntype)
            nod = nodesl(nvert(ntype)+nedge(ntype)+j)
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
#if DEBUG_MODE
            if ((RANK.eq.ROOT) .and. (iprint.eq.1)) then
               !$OMP CRITICAL
               nre = nedge(NODES(mdle)%ntype)
               nrf = nface(NODES(mdle)%ntype)
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
               !call find_element_closing_ref(ntype,kreff,krefe, kref)
!           ...-------------------------------------------------------------------
!           ...Option 2: always ask for isotropic refinement
               call get_isoref(mdle, kref)
!           ...-------------------------------------------------------------------
!           ...Option 3: always ask for radial (xy) refinement (FIBER LASER)
               !select case (NODES(mdle)%ntype)
               !   case(MDLB); kref = 110
               !   case(MDLP); kref = 10
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
#if DEBUG_MODE
      iprint = 2
      if ((RANK.eq.ROOT) .and. (iprint.eq.2)) then
         write(*,7002) ' close_mesh : number of elements to refine = ', ic
   7002  format(A,i6)
      endif
      iprint = 0
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
#if DEBUG_MODE
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

#if DEBUG_MODE
   call par_verify
#endif   
!
!
end subroutine close_mesh
