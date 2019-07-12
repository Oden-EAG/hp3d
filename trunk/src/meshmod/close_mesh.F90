subroutine close
      implicit none
      integer :: iprint
!  ...ensure backward compatibility
      iprint=0
      if (iprint /= 0)  write(*,*) 'Mesh closing begin' 
      call close_mesh
      if (iprint /= 0)  write(*,*) 'Mesh closing end'
endsubroutine close
!
!
!
!> Purpose : enforce one-irregular mesh
subroutine close_mesh()
  use error
  use refinements
  use data_structure3D
  implicit none
  ! ** Locals
  character(len=4)       :: type
7001          format('close_mesh: i= ',i6,' mdle= ', i6,' ', a4, &
                   ' ref_kind = ',i5)
  integer, allocatable   :: list(:,:)
  integer, dimension(27) :: nodesl,norientl
  integer, dimension(6)  :: kreff
  integer, dimension(12) :: krefe
  integer :: iprint, istat, i, j, ic, mdle, nod, kref, nre, nrf 
  logical :: nflag
  !-------------------------------------------------------
  iprint=0
  !
  ! allocate list
  allocate(list(2,NRELES), stat=istat)
  if (istat.ne.SUCCESS) then
    call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  endif
  !
  !
  !-------------------------------------------------------
  do 
     !---------------------------------------------------------
     ! Step 0 : activate nodes
     !---------------------------------------------------------
     call refresh
     
     !---------------------------------------------------------
     ! Step 1 : check constrained nodes
     !---------------------------------------------------------
     if (iprint.eq.2) then
        write(*,*) 'close_mesh: CHECK CONSTRAINING NODES'
     endif

     call reset_visit
     mdle=0
     do i=1,NRELES
        call nelcon(mdle, mdle)
        call get_connect_info(mdle, nodesl,norientl)
        call flag_constr_parents(mdle)
     enddo

     !---------------------------------------------------------
     ! Step 2 : pick the middle nodes to be refined
     !---------------------------------------------------------
     if (iprint.eq.2) then
        write(*,*) 'close_mesh: RESOLVE THE NODES MORE THAN ONE IRREGULARITY'
     endif

     mdle=0; ic=0
     do i=1,NRELES

        call nelcon(mdle, mdle)
        call elem_nodes(mdle, nodesl,norientl)

        type  = NODES(mdle)%type
        nflag = .FALSE.
        
        ! check edges
        !~~~~~~~~~~~~
        krefe = 0
        do j=1, nedge(type)
           nod = nodesl(nvert(type)+j)
           if (NODES(nod)%visit.eq.1) then
              krefe(j)=1
              nflag=.TRUE.
           endif
        enddo

        ! check faces
        !~~~~~~~~~~~~
        kreff = 0
        do j=1, nface(type)
           nod = nodesl(nvert(type)+nedge(type)+j)
           if (NODES(nod)%visit.eq.1) then
              call get_isoref(nod, kreff(j))
              nflag=.TRUE.
           endif
        enddo

        ! find out more than one-irregular
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (nflag) then
           if (iprint.eq.1) then
              nre = nedge(NODES(mdle)%type)
              nrf = nface(NODES(mdle)%type)
              write(*,*) 'close_mesh: mdle = ', mdle
              write(*,7003) krefe(1:nre)
7003          format('krefe = ',12i2)
              write(*,7004) kreff(1:nrf)
7004          format('kreff = ',6i2)
              call pause
           endif

           ! add counter
           ic=ic+1

           ! find out proper refinement flag
           if ( is_iso_only() ) then
              call get_isoref(mdle, kref)
           else
              call find_element_closing_ref(type,kreff,krefe, kref)
           endif
           list(1,ic) = mdle
           list(2,ic) = kref
        endif
     enddo

     ! loop out condition
     if (iprint.eq.1) then
        write(*,*) 'close_mesh: number of elements to refine ', ic
     endif
     if (ic.eq.0) then 
        exit
     end if

     !---------------------------------------------------------
     ! Step : refine from the list
     !---------------------------------------------------------
     do i=1,ic
        mdle = list(1,i)
        kref = list(2,i)
        if ( is_leaf( mdle ) ) then
           if (iprint.eq.2) then
              write(*,7001) i, mdle, NODES(mdle)%type, kref 
           endif
           call refine( mdle, kref )
        endif
     enddo
  enddo
  !
  !
  !
  deallocate(list, stat=istat) 
  if (istat.ne.SUCCESS) then
    call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  endif
!
!
endsubroutine close_mesh
