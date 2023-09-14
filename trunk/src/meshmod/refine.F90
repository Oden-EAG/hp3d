!> @brief     Refines an element enforcing TWO mesh regularity rules:
!> @details           Rule 1: no element can be broken unless
!!                            ALL its mid-face nodes are active
!!                    Rule 2: an element refinement flag is
!!                            always upgraded to accommodate
!!                            existing refinements of faces
!> @param[in] Mdle_in - middle node
!> @param[in] Kref_in - refinement kind
!> @date      Feb 2023
!--------------------------------------------------------------------
subroutine refine(Mdle_in,Kref_in)
   use data_structure3D
   use constrained_nodes
   use refinements
   use mpi_param, only: RANK,ROOT
   implicit none
!
   integer, intent(in) :: Mdle_in, Kref_in
!
   integer, dimension(MAXQUEUE) :: mdle_list,kref_list
   integer, dimension(2)        :: neig,nsid_list,norient_list
   integer, dimension(27)       :: nodesl,norientl
   integer, dimension(6)        :: kreff
   integer :: i, n, loc, loc2, iface, iface_loc, kref, krefm, kref_iso
   integer :: nrneig, nod, nodp, mdle, mdle_loc, norient_loc, nc, icase, ntype
   logical :: iflag
!
#if DEBUG_MODE
   integer :: iprint
   iprint=0
#endif
!
!---------------------------------------------------------------------
!
#if DEBUG_MODE
   if (RANK.eq.ROOT .and. iprint.ge.1) then
      write(*,*) 'refine: BEGIN Mdle_in,Kref_in,ISO = ', &
                                Mdle_in,Kref_in,is_iso_only()
   endif
#endif
!
!..if trying to refine already broken element, just return
   if (.not.is_leaf(Mdle_in)) then
      return
   endif
!
!---------------------------------------------------------------------
!  ISO ONLY - do not use shelf algorithm
!---------------------------------------------------------------------
   if ( is_iso_only() ) then
      call get_isoref(Mdle_in, krefm)
   else
      krefm = Kref_in
   endif
!
!---------------------------------------------------------------------
!  Step 0 : initialization
!---------------------------------------------------------------------
   n = 1
   mdle_list(n) = Mdle_in
   kref_list(n) = krefm
   kreff(1:6)   = 0
!
!---------------------------------------------------------------------
!  Loop over queue until queue is empty - check faces only
!---------------------------------------------------------------------
   do while (n.gt.0)
!
!  ...pick the last one
      iflag = .true.
      mdle  = mdle_list(n)
      kref  = kref_list(n)
      ntype  = NODES(mdle)%ntype
!
#if DEBUG_MODE
      if (RANK.eq.ROOT .and. iprint.ge.1) then
         write(*,*)'------------------------------------------------'
         write(*,7003) n,mdle,kref,S_Type(ntype)
         write(*,*)'------------------------------------------------'
 7003    format(' refine: n,mdle,kref,type = ',i2,',',i9,',',i3,',',a4)
      endif
#endif
!
!     Queue should not contain non-leaf element
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (.not.is_leaf(mdle)) then
         write(*,7005) mdle,NODES(mdle)%ref_kind
 7005    format('refine: ALREADY REFINED, mdle,ref_kind = ',i9,',',i3)
         call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
      endif
!
      call get_connect_info(mdle, nodesl,norientl)
!
!     loop over faces
!     ~~~~~~~~~~~~~~~~
      do iface=1,NFACE(ntype)
         i = NVERT(ntype)+NEDGE(ntype)+iface
         nod = nodesl(i)
!
#if DEBUG_MODE
         if (RANK.eq.ROOT .and. iprint.eq.2) then
            write(*,8000) iface,nod,S_Type(NODES(nod)%ntype)
8000        format('refine: active,iface,nod,type = ',i2,',',i9,',',a4)
         endif
#endif
!
!        the face is active - record the existing refinement flag
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (Is_active(nod)) then
            if (is_leaf(nod)) then
               kreff(iface) = 0
            else
               call change_ref_flag('g2l', NODES(nod)%ntype, &
                                    NODES(nod)%ref_kind,norientl(i), &
                                    kreff(iface))
            endif
!        the face is inactive (constrained)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         else
!
!        ...neighbor across the face
            call neig_face(nod, nrneig,neig,nsid_list,norient_list)
            select case (nrneig)
            case(1)
!           ...the face is boundary - record the existing refinement flag
               if (is_leaf(nod)) then
                  kreff(iface) = 0
               else
                  call change_ref_flag('g2l', NODES(nod)%ntype, &
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
!  would not work!
!  "Refesh" routine does not work with 2-irregular meshes.
!
!  Kyungjoo : refresh should be working for arbitrary irregular mesh.
!
!-----------------------------------------------------------------------
!!!              ! skip triangle
!!!              if (NODES(nod)%ntype.eq.MDLT) then
!!!                 kreff(iface) = 0
!!!                 cycle
!!!              endif
!=======================================================================
!
!           ...the face is interface of elements
               call decode2(NODES_CONSTR(i), nc,icase)
               nodp = abs(NFACEC(nc))
!
               call locate(mdle,neig,nrneig, loc)
               if (loc.eq.0) then
                  write(*,7006) mdle,iface,nod,neig
 7006             format('refine: ERROR !! mdle,iface,nod, = ',3i9, &
                        ' neig = ',2i6)
                  stop
               endif
!
!           ...get neighbor information
               loc         = mod(loc,2)+1
               mdle_loc    = neig(loc)
               iface_loc   = nsid_list(loc)
               norient_loc = norient_list(loc)
!
!           ...sanity check
               call locate(mdle_loc,mdle_list(1:n),n, loc2)
               if (loc2.gt.0) then
!              ...termination condition for the deadlock
                  write(*,*) 'refine: WARNING !! deadlock, exiting...'
                  stop
               endif
!
               kreff = 0
               call change_ref_flag('g2l',NODES(nodp)%ntype, &
                    NODES(nodp)%ref_kind,norient_loc, kreff(iface_loc))
!
!           ...-------------------------------------------------------------------
!           ...Option 1: do minimum refinement that is necessary
               !call find_element_ref(NODES(mdle_loc)%ntype,0,kreff, krefm)
!           ...-------------------------------------------------------------------
!           ...Option 2: always ask for isotropic refinement
               call get_isoref(mdle_loc, kref_iso)
               call find_element_ref(NODES(mdle_loc)%ntype,kref_iso,kreff, krefm)
!           ...-------------------------------------------------------------------
!           ...Option 3: always ask for radial (xy) refinement (FIBER LASER)
               !select case (NODES(mdle_loc)%nntype)
               !   case(MDLB); krefm = 110
               !   case(MDLP); krefm = 10
               !   case default
               !      write(*,*) 'refine: unexpected element type: mdle_loc,type = ', &
               !                          mdle_loc,S_Type(NODES(mdle_loc)%ntype)
               !end select
!           ...-------------------------------------------------------------------
!
#if DEBUG_MODE
               if (RANK.eq.ROOT .and. iprint.eq.2) then
                  write(*,8002) kreff,krefm
 8002             format(' kreff,krefm = ',6(i2,2x),',',i3)
               endif
#endif
!
               if (krefm.lt.0) then
                  write(*,*) 'refine: INCONSISTENCY'
                  write(*,*) 'mdle_loc%type,kreff,krefm      = ', &
                        S_Type(NODES(mdle_loc)%ntype),kreff,krefm
                  write(*,*) 'nodp, nodp%type, nodp%ref_kind =', &
                        nodp,S_Type(NODES(nodp)%ntype),NODES(nodp)%ref_kind
               endif
!
!           ...add neighboring element to queue
               n=n+1
               if (n.gt.MAXQUEUE) then
                  write(*,*) 'refine: MAYBE DEADLOCK '
                  call logic_error(ERR_OUT_OF_RANGE, __FILE__,__LINE__)
               endif
               mdle_list(n) = mdle_loc
               kref_list(n) = krefm
               iflag = .false.
!           ...exit loop over faces
               exit
            end select
!     ...endif node active/not active
         endif
      enddo
!
!  ...dequeue the mdle from the list
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (iflag) then
!
!         write(*,1234) 'refine: type,kref,kreff = ',S_Type(ntype),kref,kreff
! 1234    format(a,a,',',i3,',',6(i3,1x))
!
         if ( is_iso_only() ) then
            call get_isoref(mdle, krefm)
         else
            call find_element_ref(ntype,kref,kreff, krefm)
         endif
         if (kref .ne. krefm) then
            ! CHECK (FOR FIBER ONLY)
            ! fiber problem with radial refinements:
            ! krefm should always be equal to kref
         endif
!
         if (krefm.gt.0) then
#if DEBUG_MODE
            if (RANK.eq.ROOT .and. iprint.ge.1) then
               write(*,7001) mdle,kref,krefm
 7001          format('refine: BREAKING mdle,kref,krefm = ',i7,', ',i3,', ',i3)
            endif
#endif
            call break(mdle,krefm)
            n = n - 1
         else
            write(*,7010) mdle,S_Type(ntype),kref,kreff(1:NFACE(ntype)),krefm
 7010       format('refine: INCONSISTENCY, mdle,type,kref,kreff,krefm = ', &
                                           i7,',',a4,',',i3,',',6(i3,1x),',',i3)
            call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
            stop
         endif
      endif
!..end loop over queue
   enddo
!
end subroutine refine
