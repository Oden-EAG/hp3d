!> Purpose - routine refines an element enforcing TWO mesh regularity rules:
!!
!!                        Rule 1: no element can be broken unless
!!                                ALL its mid-face nodes are active
!!                        Rule 2: an element refinement flag is
!!                                always upgraded to accommodate
!!                                existing refinements of faces
!! @param[in] Mdle_in - middle node
!! @param[in] Kref_in - refinement kind
!--------------------------------------------------------------------
subroutine refine(Mdle_in,Kref_in)
  use data_structure3D
  use constrained_nodes
  use refinements
  implicit none
  ! ** Arguments
  integer, intent(in) :: Mdle_in, Kref_in
  ! ** Locals
  character(len=4) :: type
7003    format('refine: n,mdle,kref,type = ',i2,2x,i5,2x,i3,2x,a4)
 8000        format('   active,iface,nod,type = ',i2,2x,i8,2x,a4)
7010       format('refine: INCONSISTENCY, mdle,type,kref,kreff = ', &
                i5,2x,a4,2x,i2,2x,6(i2,2x))
 8001        format(' inactive,iface,nod,type = ',i2,2x,i8,2x,a4)
  !
  integer, dimension(MAXQUEUE) :: mdle_list,kref_list
  integer, dimension(2)        :: neig,nsid_list,norient_list
  integer, dimension(27)       :: nodesl,norientl
  integer, dimension(6)        :: kreff
  integer :: iprint, i, n, loc, iface, iface_loc, kref, krefm
  integer :: nrneig, nod, nodp, mdle, mdle_loc, norient_loc, nc, icase
  logical :: iflag, ideadlock
  !---------------------------------------------------------------------
  iprint=0
  if (iprint.ge.1) then
     write(*,*) 'refine: BEGIN Mdle_in,Kref_in, ISO  = ', &
          Mdle_in, Kref_in, is_iso_only()
  endif

  ! if trying to refine already broken element, just return
  if (.not.is_leaf(Mdle_in)) then
     return
  endif

  !---------------------------------------------------------------------
  ! ISO ONLY - do not use shelf algorithm
  !---------------------------------------------------------------------
  if ( is_iso_only() ) then
     call get_isoref(Mdle_in, krefm)
  else
     krefm = Kref_in
  endif

  !---------------------------------------------------------------------
  ! Step 0 : initialization
  !---------------------------------------------------------------------
  n = 1
  mdle_list(n) = Mdle_in
  kref_list(n) = krefm

  ! check it is on the shelf
  NODES(mdle_list(n))%visit = kref_list(n)

  !---------------------------------------------------------------------
  ! Loop over queue until queue is empty - check faces only
  !---------------------------------------------------------------------
  do while (n.gt.0)

     ! pick the last one
     iflag = .TRUE.
     mdle  = mdle_list(n)
     kref  = kref_list(n)
     type  = NODES(mdle)%type

     if (iprint.ge.1) then
        write(*,*)'------------------------------------------------'
        write(*,7003) n,mdle,kref,type
        write(*,*)'------------------------------------------------'
     endif

     ! Queue should not contain non-leaf element
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if (.not.is_leaf(mdle)) then
        write(*,7005) mdle,NODES(mdle)%ref_kind
7005    format('refine: ALREADY REFINED, mdle,ref_kind = ',i5,2x,i3)
        call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
     endif

     call get_connect_info(mdle, nodesl,norientl)

     ! loop over faces
     !~~~~~~~~~~~~~~~~
     do iface=1,nface(type)

        i = nvert(type)+nedge(type)+iface
        nod = nodesl(i)

        select case(NODES(nod)%act)
        case(1)

           if (iprint.eq.2) then
             write(*,8000)iface,nod,NODES(nod)%type
           endif

           ! the face is active - record the existing refinement flag
           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if (is_leaf(nod)) then
              kreff(iface) = 0
           else
              call change_ref_flag('g2l', NODES(nod)%type, &
                                   NODES(nod)%ref_kind,norientl(i), &
                                   kreff(iface))
           endif

        case(0)

            if (iprint.eq.2) then
             write(*,8000)iface,nod,NODES(nod)%type
           endif

           ! the face is inactive (constrained)
           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           ! neighbor across the face
           call neig_face(nod, nrneig,neig,nsid_list,norient_list)
           select case (nrneig)
           case(1)
              ! the face is boundary - record the existing refinement flag
              if (is_leaf(nod)) then
                 kreff(iface) = 0
              else
                 call change_ref_flag('g2l', NODES(nod)%type, &
                                      NODES(nod)%ref_kind,norientl(i), &
                                      kreff(iface))
              endif
           case(2)
!=======================================================================
!  REMARK: although there is only one refinement kind for triangular
!  faces, they cannot be skipped, otherwise
!
!    close --> refresh
!
!  would not work! "Refesh" routine does not work with 2-irregular
!  meshes.
!
!  Kyungjoo : refresh should be working for arbitrary irregular mesh.
!
!-----------------------------------------------------------------------
!!!              ! skip triangle
!!!              if (NODES(nod)%type.eq.'mdlt') then
!!!                 kreff(iface) = 0
!!!                 cycle
!!!              endif
!=======================================================================

              ! the face is interface of elements
              call decode2(NODES_CONSTR(i), nc,icase)
              nodp = abs(NFACEC(nc))

              call locate(mdle,neig,nrneig, loc)
              if (loc.eq.0) then
                 write(*,7006) mdle,iface,nod,neig
7006             format('refine: ERROR !! mdle,iface,nod, = ',3i6, &
                      ' neig = ',2i6)
                 stop 1
              endif

              ! get neighbor information
              loc         = mod(loc,2)+1
              mdle_loc    = neig(loc)
              iface_loc   = nsid_list(loc)
              norient_loc = norient_list(loc)

              if (NODES(mdle_loc)%visit.gt.0) then
                 ! ** termination condition for the deadlock
                 !    not necessary to give iso ref at this moment.
                 ! call get_isoref(mdle_list(n), kref_list(n))
                 iflag = .TRUE.
                 ! exit loop over faces
                 exit
              endif

              kreff = 0
              call change_ref_flag('g2l',NODES(nodp)%type, &
                   NODES(nodp)%ref_kind,norient_loc, kreff(iface_loc))

!!!              call find_element_ref(NODES(mdle_loc)%type,111,kreff, krefm)
              call find_element_ref(NODES(mdle_loc)%type,0,kreff, krefm)
               if (iprint.eq.2) then
                write(*,8002)kreff,krefm
 8002           format(' kreff,krefm = ',6(i2,2x),6x,2i2)
              endif

              if ((NODES(mdle_loc)%ref_filter.eq.0).and.(krefm.lt.0)) then
                 write(*,*) 'refine: INCONSISTENCY'
                 write(*,*) 'mdle_loc%type,kreff,krefm      = ', &
                      NODES(mdle_loc)%type,kreff,krefm
                 write(*,*) 'nodp, nodp%type, nodp%ref_kind =', &
                      nodp,NODES(nodp)%type,NODES(nodp)%ref_kind
                 call pause
              endif

              ! add neighboring element to queue
              n=n+1
              if (n.gt.MAXQUEUE) then
                 write(*,*) 'refine: MAYBE DEADLOCK '
                 call logic_error(ERR_OUT_OF_RANGE, __FILE__,__LINE__)
              endif
              mdle_list(n) = mdle_loc
              kref_list(n) = krefm
              NODES(mdle_list(n))%visit = kref_list(n)
              iflag = .FALSE.
              ! exit loop over faces
              exit
           end select
        end select
     enddo

     ! dequeue the mdle from the list
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if (iflag) then

        call find_element_ref(type,kref,kreff, krefm)

        if ( is_iso_only() ) then
           call get_isoref(mdle, krefm)
        endif

        if (NODES(mdle)%ref_filter.ne.0) then
           krefm = NODES(mdle)%ref_filter
        endif

        if (krefm.gt.0) then
           if (iprint.ge.1) then
              write(*,7001) mdle,kref,krefm
7001          format('refine: BREAKING mdle,kref,krefm = ',i5,2(2x,i3))
              call pause
           endif

           call break(mdle,krefm)
           NODES(mdle)%visit = 0
           n = n - 1
        else
           write(*,7010) mdle,type,kref,kreff(1:nface(type))
           call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
        endif
     endif
  enddo
!
end subroutine refine
