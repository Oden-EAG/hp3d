!-------------------------------------------------------------------------------------
!> Purpose : find YOUNGEST middle node neighbors across a mid-face node
!            who own the whole face
!!
!!
!!                   *************************************************
!!                 * *                     * *         * *         * *
!!               *   *                   *   *       *   *       *   *
!!             *     *                 *     *     *     *     *     *
!!           *************************************************       *
!!           *       *               *       *   *       *   *       *
!!           *       *               *       *************   *       *
!!           *       *               *     *\*   *     * *   *       *
!!           *       *      (1)      *   *\\\*   *   *   *   *       *
!!           *       *               * *\\\\\*   * *     *   *       *
!!           *       *               *************       *   *       *
!!           *       *               *\\\\\\\*(2)*       *   *       *
!!           *       *************************************************
!!           *     *                 *\\\\\*     *     *     *     *
!!           *   *                   *\\\*       *   *       *   *
!!           * *                     *\*         * *         * *
!!           *************************************************
!!
!!
!! @param[in]  Mface        - face node
!! @param[out] Nrneig       - number of neighbors (2 or 1)
!! @param[out] Neig         - middle node neighbors (0 if no neighbor)
!! @param[out] Nsid_list    - local face numbers in neighbors' local enumeration
!! @param[out] Norient_list - orientations of the mid-face node wrt to the neighbors
!!
!! @date Mar 2023
!-------------------------------------------------------------------------------------
subroutine neig_face(Mface, Nrneig,Neig,Nsid_list,Norient_list)
  use error
  use refinements
  use data_structure3D
!
  implicit none
! ** Arguments
  integer,               intent(in)  :: Mface
  integer,               intent(out) :: Nrneig
  integer, dimension(2), intent(out) :: Neig, Nsid_list, Norient_list
! ** Locals
  integer                    :: ntype
  integer, dimension(27, 2)  :: nodesl_neig, norientl_neig
  integer, dimension(27)     :: nodesl, norientl, nodesl_is, norientl_is
  integer, dimension(MAXGEN) :: nface_list
  integer, dimension(6)      :: nface_ort
  integer :: igen, nrgen, nve, nrf, nod, mdle, mdle_is
  integer :: kref, i, is, nrsons, iface, iflag
!
#if DEBUG_MODE
  integer :: iprint
  iprint=0
#endif
!------------------------------------------------------------------
!
!
!========================================================================
!  REMARK: 2 types of faces                                             |
!    1. faces laying on an INITIAL mesh face (they resulted from the    |
!       refinement of an inital mesh face node)                         |
!    2. faces laying inside AN element, not necessarily from the        |
!       initial mesh (they resulted from the refinement of a middle     |
!       node)                                                           |
!                                                                       |
!  If you go up the tree, you will find, for each type:                 |
!    1. the initial mesh middle node the face lays on, with mdle < 0    |
!    2. the middle node the face is inside of,         with mdle > 0    |
!       in this case, you need to go up the tree just ONCE              |
!                                                                       |
!  For each scenario:                                                   |
!    1. we use ELEMS(mdle)%neig(iface), for the proper index iface, and |
!       look for neighbors among the sons of mdle and those of its      |
!       initial mesh neighbor                                           |
!    2. we look for neighbors among the sons of mdle                    |
!========================================================================
!
  select case(NODES(Mface)%ntype)
  case(MDLT,MDLQ)
  case default
     write(*,1)Mface,S_Type(NODES(Mface)%ntype)
1    format(' neig_face: Mface,type = ',i10,',',a4)
     call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
  end select
!
  Nrneig=0 ; Neig=0 ; Nsid_list=0 ; Norient_list=0
!
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,1)Mface,S_Type(NODES(Mface)%ntype)
  endif
#endif
!
!------------------------------------------------------------------------
!  Step 1 : go UP the tree and record FACE refienement history.         |
!  REMARK: all recorded father nodes are mid-FACE nodes.                |
!------------------------------------------------------------------------
   nod=Mface ; igen=0
   do
     igen=igen+1
     nface_list(igen)=nod ; nod=NODES(nod)%father
! ...if a middle node, exit
     if (is_middle(abs(nod))) exit
  enddo
!
!  ...record number of generations
  nrgen=igen
!
#if DEBUG_MODE
  if (iprint.eq.1) then
    do i=1,nrgen
      write(*,9)i,nface_list(i)
9     format(' neig_face: i,nfather(i) = ',i2,2x,i10)
    enddo
  endif
#endif
!
  ! CASE 1 : initial mesh middle node // initial mesh face node
  ! CASE 2 : middle father node       // face node (namely Mface)
             mdle=nod                 ;  nod=nface_list(nrgen)
!
#if DEBUG_MODE
   if (iprint.eq.1) then
     write(*,3)nrgen,mdle,nod
3    format(' neig_face: nrgen,mdle,nod = ',i2,2x,2(i10,2x))
   endif
#endif
!
!..CASE 1 : initial mesh (element) middle node
   if (mdle.lt.0) then
     mdle = abs(mdle)
     ntype = ELEMS(mdle)%etype
     nve = NVERT(ntype)+NEDGE(ntype)
     nrf = NFACE(ntype)
!
! ...locate face node on list of element nodes
     call locate(nod,ELEMS(mdle)%nodes(nve+1:nve+nrf),nrf, iface)
!
! ...record "mdle" as a neighboring element
     call decodg(ELEMS(mdle)%face_orient,8,nrf, nface_ort)
     Nrneig=1 ; Neig(1)=mdle ; Nsid_list(1)=iface ; Norient_list(1)=nface_ort(iface)
!
! ...nodes and orientations of "mdle"
     call elem_nodes(mdle, nodesl_neig(1:27,1),norientl_neig(1:27,1))
!
! ...look for a neighbor of "mdle" across face
     mdle=ELEMS(mdle)%neig(iface)
!
! ...if present, record neighbor of "mdle"
     if (mdle.ne.0) then
        ntype = ELEMS(mdle)%etype
        nve = NVERT(ntype)+NEDGE(ntype)
        nrf = NFACE(ntype)
!
!    ...locate face node on list of element nodes
        call locate(nod,ELEMS(mdle)%nodes(nve+1:nve+nrf),nrf, iface)
!
        call decodg(ELEMS(mdle)%face_orient,8,nrf, nface_ort)
        Nrneig=2 ; Neig(2)=mdle ; Nsid_list(2)=iface ; Norient_list(2)=nface_ort(iface)
!
!    ...nodes and orientations of neighbor of "mdle"
        call elem_nodes(mdle, nodesl_neig(1:27,2),norientl_neig(1:27,2))
     endif
!
!..CASE 2 : internal face
   else
     ntype = NODES(mdle)%ntype
     kref = NODES(mdle)%ref_kind
     call elem_nodes(mdle, nodesl,norientl)
     call nr_mdle_sons(ntype, kref, nrsons)
     Nrneig = 0
!
! ...search for the sons
     do is=1,nrsons
        call elem_nodes_one( mdle, nodesl, norientl, is, &
                             mdle_is, nodesl_is, norientl_is )
        ntype = NODES(mdle_is)%ntype
        nve  = NVERT(ntype)+NEDGE(ntype)
        nrf  = NFACE(ntype)

        call locate(nod, nodesl_is(nve+1:nve+nrf),nrf, iface)
        if (iface.gt.0) then
           Nrneig               = Nrneig + 1
           Neig(Nrneig)         = mdle_is
           Nsid_list(Nrneig)    = iface
           Norient_list(Nrneig) = norientl_is(nve+iface)

           nodesl_neig  (1:27, Nrneig) = nodesl_is
           norientl_neig(1:27, Nrneig) = norientl_is
        endif
        if (Nrneig.eq.2) exit
     enddo
!
! ...check that two neighbors were found
     if (Nrneig.ne.2) then
       write(*,5)Mface
5      format(' neig_face: has not found 2 neighbors for internal face Mface = ',i10)
       call result
       stop
     endif
#if DEBUG_MODE
     if (iprint.eq.1) then
       write(*,4)Nrneig,Neig(1:2)
4      format(' neig_face: Nrneig,Neig(:) = ',i1,2x,2(i10,2x))
     endif
#endif
!
  endif
!
#if DEBUG_MODE
  if (iprint.eq.1) then
    write(*,7) Nrneig,Neig(1:Nrneig)
7   format(' neig_face: neighbors at top level: Nrneig,Neig = ',i1,4x,2(i10,2x))
  endif
#endif
!
!------------------------------------------------------------------------
!  Step 2: Go down the tree, sticking to the face                       |
!------------------------------------------------------------------------
!
!  igen = nrgen ; Neig = Neig(i)
!
!  DO WHILE igen > 0
!
!     face = face(igen)
!
!     look for son of Neig attached to face
!
!     IF found [update neighbor]
!
!        Neig <- son
!
!     ELSE     [advance with generation]
!
!        igen = igen - 1
!
!     ENDIF
!
!  ENDDO
!
!------------------------------------------------------------------------
!
!..loop over number of neighbors
   do i=1,Nrneig
!
     igen=nrgen
     do while (igen.ge.1)
!
        if (is_leaf(Neig(i))) then
           exit
        endif
!
!    ...pick face node
        nod =nface_list(igen)
        ntype=NODES(Neig(i))%ntype
        kref=NODES(Neig(i))%ref_kind
        call nr_mdle_sons(ntype,kref, nrsons)
#if DEBUG_MODE
        if (iprint.eq.1) then
           write(*,13)i,igen,nod,Neig(i),kref,nrsons
13         format(' neig_face: i,igen,nod,Neig,kref,nrsons = ',2(i2),2x,2i10,2i3)
           call result
        endif
#endif
!
!    ...loop over sons of Neig(i)
        iflag=0
        do is=1,nrsons
           call elem_nodes_one( Neig(i), &
                nodesl_neig(1:27,i), norientl_neig(1:27,i), is, &
                mdle_is, nodesl_is, norientl_is )
           ntype = NODES(mdle_is)%ntype
           nve = NVERT(ntype)+NEDGE(ntype)
           nrf = NFACE(ntype)
!
!       ...if face is found, update
           call locate(nod, nodesl_is(nve+1:nve+nrf),nrf, iface)
           if (iface.gt.0) then
              iflag=1
              !  ...update Neig
              Neig        (i)       = mdle_is
              Nsid_list   (i)       = iface
              Norient_list(i)       = norientl_is(nve+iface)
              nodesl_neig  (1:27,i) = nodesl_is
              norientl_neig(1:27,i) = norientl_is
!
#if DEBUG_MODE
              if (iprint.eq.1) then
                 write(*,8) i,igen,Neig(i)
8                format(' neig_face: i,igen,Neig(i) = ',i1,2x,i3,2x,i10)
              endif
#endif
              exit
           endif
!    ...end of loop over sons of Neig(i)
        enddo
        if (iflag.eq.0) igen=igen-1
     enddo
!
!  ...end of loop over Nrneig
  enddo
!
!
end subroutine neig_face
