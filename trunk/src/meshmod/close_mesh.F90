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
!
!> Purpose : enforce one-irregular mesh
subroutine close_mesh()
   use error
   use refinements
   use data_structure3D
   use mpi_param  , only: ROOT,RANK
   use MPI        , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD,MPI_REAL8
!
   implicit none
!
   character(4)           :: type
   integer                :: list(2,NRELES)
   integer                :: mdlel(NRELES)
   integer, dimension(27) :: nodesl,norientl
   integer, dimension(6)  :: kreff
   integer, dimension(12) :: krefe
   integer :: iprint, istat, i, j, ic, mdle, nod, kref, nre, nrf
   logical :: nflag
   real(8) :: MPI_Wtime,start_time,end_time
   integer :: ierr
!
#if DEBUG_MODE
   integer :: iprint
   iprint = 0
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
      if (iprint.eq.2) then
         write(*,*) 'close_mesh: CHECK CONSTRAINING NODES'
      endif
#endif
!
      list(1:2,1:NRELES) = 0; ic = 0
!
      mdle=0
      do i=1,NRELES
         call nelcon(mdle, mdle)
         mdlel(i) = mdle
      enddo
!
!$OMP PARALLEL
!$OMP DO
      do i=1,NRNODS
        NODES(i)%visit = 0
      enddo
!$OMP END DO
!$OMP DO PRIVATE(mdle,nodesl,norientl)
      do i=1,NRELES
         mdle = mdlel(i)
         call get_connect_info(mdle, nodesl,norientl)
         call flag_constr_parents(mdle)
      enddo
!$OMP END DO
!
#if DEBUG_MODE
!$OMP SINGLE
      if (iprint.eq.2) then
         write(*,*) 'close_mesh: RESOLVE THE NODES MORE THAN ONE IRREGULARITY'
      endif
!$OMP END SINGLE
#endif
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
         mdle=mdlel(i)
         call elem_nodes(mdle, nodesl,norientl)
!
         type  = NODES(mdle)%type
         nflag = .FALSE.
        
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
#if DEBUG_MODE
!$OMP CRITICAL
            if (iprint.eq.1) then
               nre = nedge(NODES(mdle)%type)
               nrf = nface(NODES(mdle)%type)
               write(*,*) 'close_mesh: mdle = ', mdle
               write(*,7003) krefe(1:nre)
7003           format('krefe = ',12i2)
               write(*,7004) kreff(1:nrf)
7004           format('kreff = ',6i2)
               call pause
            endif
!$OMP END CRITICAL
#endif
!
!           increment counter
            ic=ic+1
!
!           find out proper refinement flag
            if ( is_iso_only() ) then
               call get_isoref(mdle, kref)
            else
               call find_element_closing_ref(type,kreff,krefe, kref)
            endif
            list(1,i) = mdle
            list(2,i) = kref
         endif
      enddo
!$OMP END DO
!$OMP END PARALLEL
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,*) 'close_mesh: number of elements to refine ', ic
      endif
#endif
!
!     loop exit condition
      if (ic.eq.0) exit
!
!---------------------------------------------------------
!     Step : refine from the list
!---------------------------------------------------------
      do i=1,NRELES
         if (list(1,i) .eq. 0) cycle
         mdle = list(1,i)
         kref = list(2,i)
         if ( is_leaf( mdle ) ) then
#if DEBUG_MODE
            if (iprint.eq.2) then
               write(*,7001) i, mdle, NODES(mdle)%type, kref
 7001          format('close_mesh: i= ',i6,' mdle= ', i6,' ', a4, ' ref_kind = ',i5)
            endif
#endif
            call refine( mdle, kref )
         endif
      enddo
   enddo
!
end subroutine close_mesh
