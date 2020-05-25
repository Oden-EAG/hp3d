!-------------------------------------------------------------------------------------
!> Purpose : find YOUNGEST middle node neighbors of a mid-edge node $$$$$
!            who own the whole edge
!!
!!
!!                   *************************************************
!!                 * *                     * *         * *         * *
!!               *   *                   *   *       *   *       *   *
!!             *     *                 *     *     *     *     *     *
!!           *************************************************       *
!!           *       *               *       *   *       *   *       *
!!           *       *               *       *************   *       *
!!           *       *               *     * *   *     * *   *       *
!!           *       *      (1)      *   *   *   *   *   *   *       *
!!           *       *               * *     *   * *     *   *       *
!!           *       *               *************       *   *       *
!!           *       *               $       *(2)*       *   *       *
!!           *       ****************$********************************
!!           *     *                 $     *     *     *     *     *  
!!           *   *                   $   *       *   *       *   *    
!!           * *                     $ *         * *         * *
!!           *************************************************
!!
!!
!! @param[in]  Medge        - an edge node (marked with $$$$ above)
!! @param[in]  Maxn         - maximum number of neighbors anticipated
!! @param[out] Nrneig       - number of neighbors (>1)
!! @param[out] Neig         - middle node neighbors (0 if no neighbor)
!! @param[out] Nedg_list    - local edge numbers in neighbors' local enumeration
!! @param[out] Norient_list - orientations of the mid-edge node wrt to the neighbors
!! @param[out] Nface_list   - local face numbers in neighbors' local enumeration
!                             if the edge is contained in one of the faces
!!
!! @revision May 20
!-------------------------------------------------------------------------------------
      subroutine neig_edge(Medge,Maxn, Nrneig,Neig,Nedg_list,Norient_list,Nface_list)
!
      use error
      use refinements
      use data_structure3D
      implicit none
      common /cneig_edge/ iprint
      common /cneig_edge_mdle/ iprint_neig_edge_mdle
      common /cneig_edge_face/ iprint_neig_edge_face
      integer :: iprint_neig_edge_mdle,iprint_neig_edge_face
!
!  ...Arguments
      integer,                  intent(in)  :: Medge
      integer,                  intent(in)  :: Maxn
      integer,                  intent(out) :: Nrneig
      integer, dimension(Maxn), intent(out) :: Neig, Nedg_list, Norient_list, Nface_list
!
!  ...Locals
      character(len=4)             :: type
      integer, dimension(27,Maxn)  :: nodesl_neig, norientl_neig
      integer, dimension(27)       :: nodesl, norientl, nodesl_is, norientl_is
      integer, dimension(MAXGEN)   :: medge_ancestors
      integer, dimension(12)       :: nface_ort
      integer, dimension(Maxn)     :: neig_double
      integer :: iprint, igen, nrgen, nve, nrf, nod, mdle, mdle_is
      integer :: kref, i, is, nrsons, iface, iflag, nfath, nrv,nre, ie, loc, j
!
!========================================================================
!  REMARK: 2 types of edges                                             |
!    1. edges laying on an INITIAL mesh edge  (they resulted from the   |
!             refinement of an inital mesh edge node)                   |
!    2. edges laying inside A face, not necessarily from the            |
!       initial mesh (they resulted from the refinement of a face       |
!       node)                                                           |
!    3. edges laying inside AN element, not necessarily from the        |
!       initial mesh (they resulted from the refinement of a middle     |
!       node)                                                           |
!                                                                       |
!  If you go up the tree, you will find, for each type:                 |
!    1. the initial mesh middle node the edge lays on, with mdle < 0    |
!    2. the mid-face node the edge is inside of,       with mface > 0   |
!       in this case, you need to go up the tree just ONCE              |
!    3. the middle node the edge is inside of,         with mdle > 0    |
!       in this case, you need to go up the tree just ONCE              |
!                                                                       |
!  For each scenario:                                                   |
!    1. we use ELEMS(mdle)%neig(iface), for the proper index iface, and |
!       look for neighbors among the sons of mdle and its initial       |
!       mesh neighbors                                                  |
!    2. we look for neighbors among the sons of mface                   |
!    3. we look for neighbors among the sons of mdle                    |
!========================================================================
!
!!      select case(Medge)
!!      case(68); iprint=0
!!      case default; iprint=0
!!      end select
!
   10 continue
      iprint_neig_edge_mdle  = iprint; iprint_neig_edge_face=iprint
      select case(NODES(Medge)%type)
      case('medg')
      case default
        write(*,7010) Medge, NODES(Medge)%type
 7010   format(' neig_edge: Medge,type = ',i10,', ',a4)
        call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
      end select
!
      Nrneig=0; Neig=0; Nedg_list=0; Norient_list=0; Nface_list=0
      if (iprint.eq.1) then
        write(*,7010) Medge, NODES(Medge)%type
      endif
!
!------------------------------------------------------------------------
!  Step 1 : go UP the tree and record EDGE refinement history.         |
!  REMARK: all recorded father nodes are mid-EDGE nodes.                |     
!------------------------------------------------------------------------
      nod = Medge ; igen=0
      do
        igen=igen+1
        medge_ancestors(igen)=nod ; nod=NODES(nod)%father
!
! ......if not an edge node, exit
        if (NODES(abs(nod))%type.ne.'medg') exit
      enddo
!
!  ...record number of generations
      nrgen = igen
      if (iprint.eq.1) then
        do i=1,nrgen
          write(*,7020) i, medge_ancestors(i)
 7020     format(' neig_edge: i, medge_ancestors(i) = ',i2,2x,i10)
        enddo
      endif
!
!
      nfath = nod;  nod = medge_ancestors(nrgen) 
!
      if (iprint.eq.1) then
        write(*,7030) nrgen, nfath, nod
 7030   format(' neig_edge: nrgen, nfath, nod = ',i2,2x,2(i10,2x))
      endif
!
!  ...CASE 1: initial mesh (element) middle node
      if (nfath.lt.0) then
        call neig_edge_initial(nod,Maxn, Nrneig,Neig,nodesl_neig,norientl_neig)
      else
        select case(NODES(nfath)%type)
!
!  .....CASE 2: a face node 
        case('mdlt','mdlq')
          call neig_edge_face(nod,Maxn, Nrneig,Neig,nodesl_neig,norientl_neig,Nface_list)
!
!  .....CASE 3: a middle node 
        case('mdlp','mdlb','mdln','mdld')
          call neig_edge_mdle(nod,Maxn, Nrneig,Neig,nodesl_neig,norientl_neig)
        case default
          write(*,*) 'neig_edge: nfath, NODES(nfath)%type = ',nfath, NODES(nfath)%type 
          stop 1
        end select
      endif
!
      if (iprint.eq.1) then
        write(*,7040) nod
 7040   format(' neig_edge: NEIGHBORS FOR nod = ',i8,' AFTER STEP 1')
        do i=1,Nrneig
          type = NODES(Neig(i))%type
          nrv = nvert(type); nre = nedge(type)
          write(*,7050) Neig(i), Nodesl_neig(nrv+1:nrv+nre,i)
 7050     format(' Neig = ',i8,' Nodesl   = ',12i8)
          write(*,7060)          Norientl_neig(nrv+1:nrv+nre,i)
 7060     format('                 Norientl = ',12i8)
        enddo
      endif
!
!  ...determine local edge number
      do i=1,Nrneig
        type=NODES(Neig(i))%type
        nrv =nvert(type); nre = nedge(type)
        call locate(nod, nodesl_neig(nrv+1:nrv+nre,i),nre, ie)
        Nedg_list(i) = ie
        if (ie.ne.0) Nface_list(i)=0
      enddo
!
!------------------------------------------------------------------------
!  Step 2: Go down the tree, sticking to the edge                       |
!------------------------------------------------------------------------
!
!  ...check for double neighbors (may originate from a face ancestor)
      neig_double=1
      do i=2,Nrneig
        call locate(Neig(i), Neig(1:i-1),i-1, loc)
        if (loc.ne.0) then
          neig_double(i)=2
!
!  .......check consistency
          if ((Nedg_list(loc).ne.0).or.(Nedg_list(i).ne.0)) then
            write(*,*) 'neig_edge: INCONSISTENCY 1'
            iprint=1
            go to 10
          endif
        endif
      enddo
      if (iprint.eq.1) then
        write(*,7065) neig_double(1:Nrneig)
 7065   format(' neig_edge: neig_double = ',20i2)
      endif
!
!  ...loop over number of neighbors
      do i=1,Nrneig
!
        igen = nrgen
        do while (igen.ge.1)
!
          if (is_leaf(Neig(i))) then
            exit
          endif
!
!  .......pick edge node
          nod  = medge_ancestors(igen)
          type = NODES(Neig(i))%type
          kref = NODES(Neig(i))%ref_kind
          call nr_mdle_sons(type, kref, nrsons)
!        
          if (iprint.eq.1) then
            write(*,7070) i,igen,nod,Neig(i),kref,nrsons
 7070       format(' neig_edge: i,igen =', 2i2,' nod,Neig = ',2i6,' kref,nrsons = ',i3,i2)
            write(*,7075) Neig(i)
 7075       format(' Neig(i) = ',i6,' nodesl,norientl = ')
            write(*,7076) nodesl_neig(:,i)
 7076       format(27i5)
            write(*,7076) norientl_neig(:,i)
          endif
!
!  .......loop over sons of Neig(i)
          iflag=0
          do is=1,nrsons
            call elem_nodes_one( Neig(i), &
                                 nodesl_neig(:,i), norientl_neig(:,i), is, &
                                 mdle_is, nodesl_is, norientl_is )
            type=NODES(mdle_is)%type
            nrv =nvert(type); nre = nedge(type)
!
!  .........if edge is found, update
            call locate(nod, nodesl_is(nrv+1:nrv+nre),nre, ie)
            if (ie.gt.0) then
              iflag = iflag+1
              if (iflag.eq.neig_double(i)) then
!
!  .............a one time action, eliminate a repetition
                neig_double(i)=1
!
!  .............update Neig
                Neig(i)            = mdle_is
                Nedg_list(i)       = ie
                Norient_list(i)    = norientl_is(nrv+ie)
                nodesl_neig(:,i)   = nodesl_is
                norientl_neig(:,i) = norientl_is
!
                if (iprint.eq.1) then
                  write(*,7080) i,igen,Neig(i)
 7080             format(' neig_edge: i,igen,Neig(i) = ',i2,2x,i3,2x,i10)
                endif
                exit
              endif
            endif
!
!  .......end of loop over sons of Neig(i)
          enddo
          if (iflag.eq.0) igen=igen-1
        enddo
!
!  ...end of loop over Nrneig
      enddo
!
!  ...eliminate possible duplication
      i=2
      do while (i.le.Nrneig)
        call locate(Neig(i), Neig(1:i-1),i-1, loc)
        if (loc.ne.0) then
!
!  .......check consistency
          if ((Nedg_list(loc).ne.0).or.(Nface_list(loc).eq.0)) then
            write(*,*) 'neig_edge: INCONSISTENCY 2'
            iprint=1
            write(*,7090) Medge
            do i=1,Nrneig
              write(*,7100) Neig(i),Nedg_list(i),Norient_list(i),Nface_list(i)
            enddo
            go to 10
          endif
!
!  .......restact the output to avoid the repetition
          Nrneig = Nrneig-1
          do j=i,Nrneig
            Neig(j) = Neig(j+1)
            Nedg_list(j) = Nedg_list(j+1)
            Norient_list(j) = Norient_list(j+1)
            Nface_list(j) =  Nface_list(j+1)
          enddo
        endif
        i=i+1
      enddo
!
!
      if (iprint.eq.1) then
        write(*,7090) Medge
 7090   format(' neig_edge: Neig,Nedg,Norient,Nface FOR Medge = ',i6)
        do i=1,Nrneig
          write(*,7100) Neig(i),Nedg_list(i),Norient_list(i),Nface_list(i)
 7100     format(i8,i3,2i3)
        enddo
        call pause
      endif
!
!
      end subroutine neig_edge


!-------------------------------------------------------------------------------------
!> Purpose : find middle node neighbors of an initial mesh mid-edge node
!!
!! @param[in]  Medge        - an edge node 
!! @param[in]  Maxn         - maximum number of neighbors anticipated
!! @param[out] Nrneig       - number of neighbors (>1)
!! @param[out] Neig         - middle node neighbors (0 if no neighbor)
!! @param[out] Nodesl_neig   - nodes for the neighbors
!! @param[out] Norientl_neig - node orientations for the neighbors
!!
!! @revision May 20
!-------------------------------------------------------------------------------------
      subroutine neig_edge_initial(Medge,Maxn, Nrneig,Neig,Nodesl_neig,Norientl_neig)
!
      use GMP
      use data_structure3D
      use element_data
      implicit none
!
!  ...Arguments
      integer,                     intent(in)  :: Medge
      integer,                     intent(in)  :: Maxn
      integer,                     intent(out) :: Nrneig
      integer, dimension(Maxn),    intent(out) :: Neig
      integer, dimension(27,Maxn), intent(out) :: Nodesl_neig,Norientl_neig
!
!  ...neighbors of a curve
      integer, dimension(Maxn)  ::neigbl
!
!  ...element/node type
      character(len=4) :: type
!
      integer :: iprint, nod, nc, nrbl, ib, nb, lab, nel, nrn, nrv, nre, i
!
!----------------------------------------------------------------------
!
      iprint=0
!
      if (NODES(Medge)%type.ne.'medg') then
        write(*,7010) Medge, NODES(Medge)%type
 7010   format('neig_edge_initial: Medge, NODES(Medge)%type = ',i6,2x,a5)
        stop 1
      endif
!
!  ...consistency check
      if (NODES(Medge)%father.gt.0) then
        write(*,7020) Medge,NODES(Medge)%father
 7020   format('neig_edge_initial: Medge, NODES(Medge)%father = ',i6,2x,i6)
        stop 1
      endif
!
!  ...the GMP curve number is (see hp3gen)
      nod = Medge
      nc = nod-(NRELIS+NRPOINT)
!
!  ...use GMP utility to find the adjacent GMP blocks
      call find_curve_to_block(nc,Maxn, nrbl,neigbl)
      if (nrbl.gt.Maxn) then
        write(*,7030) Maxn,nrbl
 7030   format('neig_edge_initial: Maxn,nrbl = ',2i4)
        stop 1
      endif
      Nrneig = nrbl
!
!  ...loop through neighboring blocks
      do ib=1,nrbl
        call decode(neigbl(ib), nb,lab)
!
!  .....determine element number
        select case(lab)
!
!  .....prism
        case(1); nel= nb
!
!  .....hexahedron
        case(2); nel= NRPRISM+nb
!
!  .....tetrahedron
        case(3); nel= NRPRISM+NRHEXAS+nb
!
!  .....pyramid
        case(4); nel= NRPRISM+NRHEXAS+NRTETRA+nb
        end select
!
!  .....store the neighbor, for the initial mesh
!       middle node numbers coincide with element numbers...
        Neig(ib) = nel
        call elem_nodes(nel, Nodesl_neig(:,ib),Norientl_neig(:,ib))
!
!  ...end of loop through adjacent blocks
      enddo
!
      if (iprint.eq.1) then
        write(*,7050) Medge
 7050   format(' neig_edge_initial: NEIGHBORS FOR Medge = ',i8)
        do i=1,Nrneig
          type = NODES(Neig(i))%type
          nrv = nvert(type); nre = nedge(type)
          write(*,7060) Neig(i), Nodesl_neig(nrv+1:nrv+nre,i)
 7060     format(' Neig = ',i8,' Nodesl   = ',12i8)
          write(*,7070)          Norientl_neig(nrv+1:nrv+nre,i)
 7070     format('                 Norientl = ',12i8)
        enddo
        call pause
      endif
!
!
      end subroutine neig_edge_initial



!-------------------------------------------------------------------------------------
!> Purpose : find middle node neighbors for an edge son of a face node
!!
!! @param[in]  Medge         - an edge node 
!! @param[in]  Maxn          - maximum number of neighbors anticipated
!! @param[out] Nrneig        - number of neighbors (>1)
!! @param[out] Neig          - middle node neighbors
!! @param[out] Nodesl_neig   - nodes for the neighbors
!! @param[out] Norientl_neig - node orientations for the neighbors
!! @param[out] Nface_list   - local face numbers in neighbors' local enumeration
!                             if the edge is contained in one of the faces
!!
!! @revision May 20
!-------------------------------------------------------------------------------------
      subroutine neig_edge_face(Medge,Maxn, Nrneig,Neig,Nodesl_neig,Norientl_neig,Nface_neig)
!
      use data_structure3D
      use element_data
      implicit none
      common /cneig_edge_face/ iprint
      common /cneig_face_extended/ iprint_neig_face_extended
      integer :: iprint_neig_face_extended
!
!  ...Arguments
      integer,                     intent(in)  :: Medge
      integer,                     intent(in)  :: Maxn
      integer,                     intent(out) :: Nrneig
      integer, dimension(Maxn),    intent(out) :: Neig,Nface_neig
      integer, dimension(27,Maxn), intent(out) :: Nodesl_neig,Norientl_neig
!
!  ...work space for neig_face_extended
      integer                :: mface
      integer                :: nrneig_face
!
!  ...element/node type
      character(len=4) :: type
!
!  ...brother nodes
      integer, dimension(2)  :: brother_mface
!
      integer :: iprint, nfath, is, i, nrv, nre
!
!----------------------------------------------------------------------
!
!!!      select case(Medge)
!!!      case(68); iprint=0
!!!      case default; iprint=0
!!!      end select
!
      iprint_neig_face_extended = iprint
      if (NODES(Medge)%type.ne.'medg') then
        write(*,7010) Medge, NODES(Medge)%type
 7010   format(' neig_edge_face: Medge, NODES(Medge)%type = ',i6,2x,a5)
        stop 1
      endif
      if (iprint.eq.1) then
        write(*,7015) Medge
 7015   format(' neig_edge_face: DEBUGGING FOR Medge = ',i6)
      endif
!
!  ...father node
      nfath = NODES(Medge)%father
!
!  ...determine the son number
      do is=1,NODES(nfath)%nr_sons
        if (son(nfath,is).eq.Medge) exit
      enddo
!
!  ...determine brother face nodes
      select case(NODES(nfath)%type)
!
      case('mdlt')
        select case(NODES(nfath)%ref_kind)
        case(4)
          brother_mface(1)=is-4; brother_mface(2)=4   
        case default
          write(*,*) 'neig_edge_face: UNSUPPORTED 1'
          stop 1
        end select
      case('mdlq')
        select case(NODES(nfath)%ref_kind)
        case(11)
          brother_mface(1)=is-4; brother_mface(2)=is-3
          if (brother_mface(2).eq.5)  brother_mface(2)=1
        case(10,01)  
          brother_mface(1)=1; brother_mface(2)=2
        case default
          write(*,*) 'neig_edge_face: UNSUPPORTED 2'
          stop 1
        end select
      case default
        write(*,*) 'neig_edge_face: NODES(nfath)%type = ',NODES(nfath)%type
        stop 1
      end select
      if (iprint.eq.1) then
        write(*,7016) is, brother_mface(1:2)
 7016   format(' neig_edge_face: is = ',i2,' brother_mface = ',2i3)
      endif
!
!  ...loop through the brother face nodes
      Nrneig=0
      do i=1,2
        mface = son(nfath,brother_mface(i))
!
!  .....determine the corresponding middle node neighbors
        call neig_face_extended(mface, nrneig_face,Neig(Nrneig+1:Nrneig+2),Nface_neig(Nrneig+1:Nrneig+2), &
                                Nodesl_neig(:,Nrneig+1:Nrneig+2),Norientl_neig(:,Nrneig+1:Nrneig+2))
        Nrneig = Nrneig + nrneig_face
      enddo
!
      if (iprint.eq.1) then
        write(*,7020) Medge
 7020   format(' neig_edge_face: NEIGHBORS FOR Medge = ',i8)
        do i=1,Nrneig
          type = NODES(Neig(i))%type
          nrv = nvert(type); nre = nedge(type)
          write(*,7030) Neig(i), Nodesl_neig(nrv+1:nrv+nre,i)
 7030     format(' Neig = ',i8,' Nodesl   = ',12i8)
          write(*,7040)          Norientl_neig(nrv+1:nrv+nre,i)
 7040     format('                 Norientl = ',12i8)
        enddo
        call pause
      endif
!
!
      end subroutine neig_edge_face


          
!-------------------------------------------------------------------------------------
!> Purpose : find middle node neighbors for an edge son of a middle node
!!
!! @param[in]  Medge         - an edge node 
!! @param[in]  Maxn          - maximum number of neighbors anticipated
!! @param[out] Nrneig        - number of neighbors 
!! @param[out] Neig          - middle node neighbors
!! @param[out] Nodesl_neig   - nodes for the neighbors
!! @param[out] Norientl_neig - node orientations for the neighbors
!!
!! @revision May 20
!-------------------------------------------------------------------------------------
      subroutine neig_edge_mdle(Medge,Maxn, Nrneig,Neig,Nodesl_neig,Norientl_neig)
!
      use data_structure3D
      use element_data
      use refinements
      implicit none
      common /cneig_edge_mdle/ iprint
!
!  ...Arguments
      integer,                     intent(in)  :: Medge
      integer,                     intent(in)  :: Maxn
      integer,                     intent(out) :: Nrneig
      integer, dimension(Maxn),    intent(out) :: Neig
      integer, dimension(27,Maxn), intent(out) :: Nodesl_neig,Norientl_neig
!
!  ...element/node type
      character(len=4) :: type
!
!  ...nodes and their orientation 
      integer, dimension(27)  :: nodesl,norientl,nodesl_is,norientl_is
!
      integer :: iprint, nfath, is, i, nrsons, mdle_is, nrv, nre, loc
!
!----------------------------------------------------------------------
!
!!!      iprint=0
!
      if (NODES(Medge)%type.ne.'medg') then
        write(*,7010) Medge, NODES(Medge)%type
 7010   format(' neig_edge_mdle: Medge, NODES(Medge)%type = ',i6,2x,a5)
        stop 1
      endif
!
!  ...father node
      nfath = NODES(Medge)%father
!
!  ...determine nodes of the father and their orientation
      call elem_nodes(nfath, nodesl, norientl)
!
!  ...determine the number of middle node sons of the father
      call nr_mdle_sons(NODES(nfath)%type,NODES(nfath)%ref_kind, nrsons)
!
!  ...loop through the middle node sons
      Nrneig = 0
      do is=1,nrsons
!
!  .....determine nodes of the sons 
        call elem_nodes_one(nfath, nodesl, norientl, is, &
                            mdle_is, nodesl_is, norientl_is)
        type = NODES(mdle_is)%type
        nrv  = nvert(type); nre = nedge(type)
        call locate(Medge, nodesl_is(nrv+1:nrv+nre),nre, loc)
        if (loc.gt.0) then
           Nrneig               = Nrneig + 1
           Neig(Nrneig)         = mdle_is
           Nodesl_neig(:,Nrneig)   = nodesl_is
           Norientl_neig(:,Nrneig) = norientl_is
        endif
      enddo
!
      if (iprint.eq.1) then
        write(*,7020) Medge
 7020   format(' neig_edge_mdle: NEIGHBORS FOR Medge = ',i8)
        do i=1,Nrneig
          type = NODES(Neig(i))%type
          nrv = nvert(type); nre = nedge(type)
          write(*,7030) Neig(i), Nodesl_neig(nrv+1:nrv+nre,i)
 7030     format(' Neig = ',i8,' Nodesl   = ',12i8)
          write(*,7040)          Norientl_neig(nrv+1:nrv+nre,i)
 7040     format('                 Norientl = ',12i8)
        enddo
        call pause
      endif
!
!
      end subroutine neig_edge_mdle



          
