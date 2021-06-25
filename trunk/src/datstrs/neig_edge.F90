!-------------------------------------------------------------------------------------
!> Purpose : find YOUNGEST middle node neighbors along a mid-edge node
!            who own the whole edge
!!
!!
!! @param[in]  Medge        - an edge node
!! @param[in]  Maxn         - maximum number of neighbors 
!! @param[out] Nrneig       - number of neighbors 
!! @param[out] Neig         - middle node neighbors
!! @param[out] Nedg_list    - local edge numbers in neighbors' local enumeration
!! @param[out] Norient_list - orientations of the mid-edge node wrt to the neighbors
!!
!! @revision Jun 21
!-------------------------------------------------------------------------------------
!
      subroutine neig_edge(Medge,Maxn, Nrneig,Neig,Nedg_list,Norient_list)
      use error
      use data_structure3D
!
      implicit none
!
!  ...Arguments
      integer, intent(in)  :: Medge
      integer, intent(in)  :: Maxn
      integer, intent(out) :: Nrneig
      integer, intent(out) :: Neig(Maxn), Nedg_list(Maxn), Norient_list(Maxn)
!
!  ...Locals
      integer :: nedg_ancestors(MAXGEN)
      integer :: iprint,nod,igen,nrgen,nfath,nson,i,j,nrv,nre,loc  
!
!  ...element nodes and orientation
      integer :: nodesl_fath(27),norientl_fath(27)
      integer :: nodesl_son(27), norientl_son(27)
      integer :: neig_nodesl(27,Maxn),neig_orientl(27,Maxn)
      integer :: nr_mdle_sons
!
!
!========================================================================
!  REMARK: 3 types of edges:                                            |
!    1. edges laying on an INITIAL mesh edge, they result from          |
!       refinements of the initial mesh edge                            |
!    2. edges laying on an INITIAL mesh face (they result from          |
!       refinements of the inital mesh face node)                       |
!    2. faces laying inside AN element, not necessarily from the        |
!       initial mesh (they resulted from the refinement of a middle     |
!       node)                                                           |
!========================================================================
!
      select case(Medge)
      case(139)
        iprint=0
      case default
        iprint=0
      end select
!
      select case(NODES(Medge)%type)
      case('medg')
      case default
        write(*,7100) Medge,NODES(Medge)%type
 7100   format(' neig_edge: Medge,type = ',i10,',',a4)
        call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
      end select
!
      Nrneig=0 ; Neig=0 ; Nedg_list=0 ; Norient_list=0
      if (iprint.eq.1) then
        write(*,7100) Medge,NODES(Medge)%type
      endif
!
!------------------------------------------------------------------------
!  Step 1 : go UP the tree and record edge refinement history.         |
!------------------------------------------------------------------------
      nod=Medge ; igen=0
      do
        igen=igen+1
        nedg_ancestors(igen)=nod 
        nfath = NODES(nod)%father
!
!  .....we have reached an initial mesh element edge
        if (nfath.lt.0) then
!
!  .......determine initial mesh element neighbors of the edge
          call neig_initial_mesh_edge(nod,Maxn, Nrneig, &
                                      Neig,neig_nodesl,neig_orientl)
          go to 200
!
        else
          select case(NODES(nfath)%type)
!
!  .......father is an edge, proceed up the tree
          case('medg')
            nod = nfath
            go to 100
!
!  .......father is a mid-face node
          case('mdlt','mdlq')
            call neig_mface_edge(nod,Maxn, Nrneig, &
                                 Neig,neig_nodesl,neig_orientl)
            go to 200
!
!  .......father is a middle node, determine the edge neighbors
!         from the list of its sons
          case('mdlb','mdln','mdlp','mdld')
            call neig_middle_edge(nod,Maxn, Nrneig, &
                                  Neig,neig_nodesl,neig_orientl)
            go to 200
          end select

  100   continue
        endif
      enddo
!
!  ...Step 2: Go down the tree...
  200 continue
      nrgen = igen
      if (iprint.eq.1) then
        write(*,7030) nrgen
 7030   format('neig_edge: nrgen = ',i3)
        write(*,7040) Neig(1:Nrneig)
 7040   format('neig_edge: NEIGBORS AFTER THE FIRST STEP = ',10i10)
        call pause
      endif
!
!  ...loop through the neighbors of the ancestor edge
      do i=1,Nrneig
        nfath = Neig(i); igen = nrgen 
        nod = nedg_ancestors(igen) 
!
!  .....recover nodes for the middle node element
        nodesl_fath(1:27) = neig_nodesl(1:27,i)
        norientl_fath(1:27) = neig_orientl(1:27,i)
  300   Neig(i) = nfath
        if (NODES(nfath)%ref_kind.eq.0) then
          if (igen.ne.1) then
            write(*,*) 'neig_edge: INCONSISTENCY, igen = ',igen
            stop 1
          endif
!
!  .......look for the current nod in the list of the element nodes
          nrv = nvert(NODES(nfath)%type)
          nre = nedge(NODES(nfath)%type)
          call locate(nod,nodesl_fath(nrv+1:nrv+nre),nre, loc)
          go to 400
        endif
!
!  .....loop through the middle node sons of the mdle node
        do j=1,nr_mdle_sons(NODES(nfath)%type,NODES(nfath)%ref_kind)
!
!  .......determine nodal connectivities for the son
          call elem_nodes_one(nfath,nodesl_fath,norientl_fath,j, &
                              nson, nodesl_son, norientl_son )
!
!  .......look for the current nod in the list of the element nodes
          nrv = nvert(NODES(nson)%type)
          nre = nedge(NODES(nson)%type)
          call locate(nod,nodesl_son(nrv+1:nrv+nre),nre, loc)
!
!  .......if you have found the current node on the list, replace
!         father with the son
          if (loc.gt.0) then
            nfath = nson
            nodesl_fath = nodesl_son
            norientl_fath = norientl_son
            go to 300
          endif
!
!  .......look for the son of the current edge in the list of the element nodes
          call locate(nedg_ancestors(igen-1),nodesl_son(nrv+1:nrv+nre),nre, loc)
          if (iprint.eq.1) then
            write(*,*) 'i,j,loc = ',i,j,loc
          endif
!
!  .......if you have found the node on the list, replace
!         father with the son
          if (loc.gt.0) then
            nod = nedg_ancestors(igen-1)
            nfath = nson
            nodesl_fath = nodesl_son
            norientl_fath = norientl_son
!
!  .........update the generation number and switch to the son of 'nod'
            igen = igen-1
            go to 300
          endif
        enddo
  400   continue
!
!  .....update the edge number and orientation of the node
        Nedg_list(i) = loc
        Norient_list(i) = norientl_fath(nrv+loc)
!
!  ...end of loop through the neighbors
      enddo
!
      if (iprint.eq.1) then
        write(*,7110) Medge
 7110   format('neig_edge: Neig,Nedg_list,Norient_list FOR Medge = ',i10)
        write(*,7020) Neig(1:Nrneig)
 7020   format(8i10)
        write(*,7020) Nedg_list(1:Nrneig)
        write(*,7020) Norient_list(1:Nrneig)
        call pause
      endif
!
      end subroutine neig_edge


!----------------------------------------------------------------------
!
!   routine name       - neig_initial_mesh_edge
!
!----------------------------------------------------------------------
!
!   latest revision    - Jun 21
!
!   purpose            - find neigbors for an INITIAL mesh element edge
!
!   arguments :
!     in:
!         Medg         - a mid-edge node
!         Maxn         - dimension of the arrays below
!                        (max number of neighbors)
!     out:
!         Nrneig       - number of adjacent elements
!         Neig         - list of neighbors of the edge
!         Neig_nodesl,Neig_orientl - nodal connectivities for the
!                        neighbors
!
!   required  routines -
!
!----------------------------------------------------------------------
!   
      subroutine neig_initial_mesh_edge(Medg,Maxn, Nrneig, &
                                        Neig,Neig_nodesl,Neig_orientl)
!
      use GMP
      use data_structure3D
      use element_data
      implicit none
!
!  ...Arguments
      integer, intent(in)  :: Medg
      integer, intent(in)  :: Maxn
      integer, intent(out) :: Nrneig
      integer, intent(out) :: Neig(Maxn)
      integer, intent(out) :: Neig_nodesl(27,Maxn),Neig_orientl(27,Maxn)
!
!  ...neighbors of a curve
      integer :: neigbl(100)
!
!  ...decoded orientations for edges
      integer :: nedge_orient(12)
!
!  ...element/node type
      character(len=4) :: type
!
      integer :: nod,nc,nrbl,ib,nb,lab,nel,iprint
!
!----------------------------------------------------------------------
!
      iprint=0
!
      if (NODES(Medg)%type.ne.'medg') then
        write(*,7001) Medg,NODES(Medg)%type
 7001   format('neig_initial_mesh_edge: Medg,NODES(Medg)%type = ',i6,2x,a5)
        stop 1
      endif
!
!  ...a version for the initial mesh only...
      if (NODES(Medg)%father.gt.0) then
        write(*,7002) Medg,NODES(Medg)%father
 7002   format('neig_initial_mesh_edge: Medg,NODES(Medg)%father = ',i6,2x,i6)
        stop 1
      endif
!
!  ...the number of the GMP curve is (see hp3gen)
      nod = Medg
      nc = nod-(NRELIS+NRPOINT)
!
!  ...use GMP utility to find the adjacent GMP blocks
      call find_curve_to_block(nc,100, nrbl,neigbl)
      if (nrbl.gt.Maxn) then
        write(*,7003) Maxn,nrbl
 7003   format('neig_initial_mesh_edge: Maxn,nrbl = ',2i4)
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
        call elem_nodes(nel,Neig_nodesl(1:27,ib),Neig_orientl(1:27,ib))
!
!  ...end of loop through adjacent blocks
      enddo
      if (iprint.eq.1) then
        write(*,7011) Medg
 7011   format('neig_initial_mesh_edge: NEIGHBORS OF Medg = ',i5)
        write(*,7012) Neig(1:Nrneig)
 7012   format(10i8)
        call pause
      endif
!
!
      end subroutine neig_initial_mesh_edge

!----------------------------------------------------------------------
!
!   routine name       - neig_middle_edge
!
!----------------------------------------------------------------------
!
!   latest revision    - Jun 21
!
!   purpose            - find midle node neigbors for an edge son
!                        of a middle node
!
!   arguments :
!     in:
!         Medg         - a mid-edge node
!         Maxn         - dimension of the arrays below
!                        (max number of neighbors)
!     out:
!         Nrneig       - number of adjacent elements
!         Neig         - list of neighbors of the edge
!         Nedg_list    - local edge numbers
!         Norient_list - orientations of the mid-edge node
!
!----------------------------------------------------------------------
!   
      subroutine neig_middle_edge(Medg,Maxn, Nrneig, &
                                  Neig,Neig_nodesl,Neig_orientl)
!
      use GMP
      use data_structure3D
      use element_data
      implicit none
!
!  ...Arguments
      integer, intent(in)  :: Medg
      integer, intent(in)  :: Maxn
      integer, intent(out) :: Nrneig
      integer, intent(out) :: Neig(Maxn)
      integer, intent(out) :: Neig_nodesl(27,Maxn),Neig_orientl(27,Maxn)
!
!  ...element nodes and orientation
      integer :: nodesl_fath(27),norientl_fath(27)
      integer :: nodesl_son(27), norientl_son(27)
!
      integer :: nr_mdle_sons
      integer :: mdle,nrsons,is,nod,nrv,nre,loc,i,iprint,nson
!
!----------------------------------------------------------------------
!
      iprint=0
!
      mdle = NODES(Medg)%father
      select case(NODES(mdle)%type)
      case('mdlb','mdln','mdlp','mdld')
      case default
        write(*,7100) Medg
 7100   format('neig_middle_edge: WRONG FATHER OF Medg = ',i10)
        stop 1
      end select 
      Neig = 0
!
!  ...determine nodes for the middle node element
      call elem_nodes(mdle, nodesl_fath,norientl_fath) 
!
!  ...initiate number of neighbors
      i=0
!
!  ...loop through the middle node sons
      do is=1,nr_mdle_sons(NODES(mdle)%type,NODES(mdle)%ref_kind)
        nson = Son(Mdle,is)
!
!  .....determine nodes for the son using the nodes of the father
        call elem_nodes_one(mdle,nodesl_fath,norientl_fath,is, &
                            nod,nodesl_son,norientl_son)
        nrv = nvert(NODES(nod)%type)
        nre = nedge(NODES(nod)%type)
        call locate(Medg,nodesl_son(nrv+1:nrv+nre),nre, loc)
        if (loc.gt.0) then
          i=i+1
          if (i.gt.Maxn) then
            write(*,*) 'neig_middle_edge: INSUFFIOCIENT Maxn = ',Maxn
            stop 1
          endif
          Neig(i) = nod
          Neig_nodesl(1:27,i)  = nodesl_son(1:27)
          Neig_orientl(1:27,i) = norientl_son(1:27)
        endif
!
!  ...end of loop through the middle node sons
      enddo
      Nrneig = i
!

      if (iprint.eq.1) then
        write(*,7110) Medg
 7110   format('neig_middle_edge: NEIGHBORS OF Medg = ',i5)
        write(*,7120) Neig(1:Nrneig)
 7120   format(10i8)
        call pause
      endif
!
!
      end subroutine neig_middle_edge


!----------------------------------------------------------------------
!
!   routine name       - neig_mface_edge
!
!----------------------------------------------------------------------
!
!   latest revision    - Jun 21
!
!   purpose            - find midle node neigbors for an edge son
!                        of a face node
!
!   arguments :
!     in:
!         Medg         - a mid-edge node
!         Maxn         - dimension of the array below
!                        (max number of neighbors)
!     out:
!         Nrneig       - number of adjacent elements
!         Neig         - list of neighbors of the edge
!         Neig_nodesl,Neig_orientl - nodal connectivities for the
!                        neighbors
!
!----------------------------------------------------------------------
!   
      subroutine neig_mface_edge(Medg,Maxn, Nrneig,&
                                 Neig,Neig_nodesl,Neig_orientl)
!
      use data_structure3D
      use element_data
      implicit none
!
!  ...Arguments
      integer, intent(in)  :: Medg
      integer, intent(in)  :: Maxn
      integer, intent(out) :: Nrneig
      integer, intent(out) :: Neig(Maxn)
      integer, intent(out) :: Neig_nodesl(27,Maxn),Neig_orientl(27,Maxn)
!
      integer :: mface,is,medg_bro(2),i,j,loc,nrv,nre
      integer, parameter, dimension(1:2,1:3) :: mdlt4_bro  = reshape((/1,4, 2,4, 3,4/), (/2,3/))
      integer, parameter, dimension(1:2,1:4) :: mdlq11_bro = reshape((/1,2, 2,3, 3,4, 4,1/), (/2,4/))
!
!  ...work space for routine 'neig_face'
      integer :: nrneigf,neigf(2),nsidf(2),norientf(2)
!
!  ...element nodes and orientation
      integer :: nodesl(27),norientl(27)
!
      integer :: iprint
!
!----------------------------------------------------------------------
!
      iprint=0
!
      mface = NODES(Medg)%father
      select case(NODES(mface)%type)
      case('mdlt','mdlq')
      case default
        write(*,7100) Medg
 7100   format('neig_mface_edge: WRONG FATHER OF Medg = ',i10)
        stop 1
      end select 
      Nrneig=0
      Neig = 0
!
!  ...locate 'Medg' on the list of nodal sons of 'Mface'
      do is = 1,NODES(mface)%nr_sons
        if (Medg.eq.son(mface,is)) exit
      enddo
!
!  ...determine mid-face brothers of 'Medg'
      do i=1,2
        select case(NODES(mface)%type)
        case('mdlt')
          medg_bro(i) = son(mface,mdlt4_bro(i,is-4))
        case('mdlq')
          select case(NODES(mface)%ref_kind)
          case(11)
            medg_bro(i) = son(mface,mdlq11_bro(i,is-4)) 
          case(10,01)
            medg_bro(i) = son(mface,i) 
          end select
        end select
      enddo
!
!  ...determine face neighbors for each of face brothers
      do i=1,2
        call neig_face(medg_bro(i), nrneigf,neigf,nsidf,norientf)
        do j=1,nrneigf
!
!  .......determine nodes of the neighbor
          call elem_nodes(neigf(j), nodesl,norientl)
          nrv = nvert(NODES(neigf(j))%type)
          nre = nedge(NODES(neigf(j))%type)
          call locate(Medg, nodesl(nrv+1:nrv+nre),nre, loc)
          if (loc.gt.0) then
!
!  .........add the element middle node to the list
            Nrneig = Nrneig+1
            if (Nrneig.gt.Maxn) then
              write(*,*) 'neig_mface_edge: INSUFFIOCIENT Maxn = ',Maxn
              stop 1
            endif
            Neig(Nrneig) = neigf(j)
            Neig_nodesl(1:27,Nrneig)  = nodesl(1:27)
            Neig_orientl(1:27,Nrneig) = norientl(1:27)
          endif
        enddo
      enddo
!
      if (iprint.eq.1) then
        write(*,7110) Medg
 7110   format('neig_mface_edge: NEIGHBORS OF Medg = ',i5)
        write(*,7120) Neig(1:Nrneig)
 7120   format(10i8)
        call pause
      endif
!
!
      end subroutine neig_mface_edge


!----------------------------------------------------------------------
!
!   function name      - Nr_mdle_sons
!
!----------------------------------------------------------------------
!
!   latest revision    - Jun 21
!
!   purpose            - return number of middle node sons for a middle
!                        node
!
!   arguments :
!     in:
!         Type         - middle node type
!         Kref         - refinement flag
!     out:
!         Nr_mdle_sons - number of middle nod sons
!
!----------------------------------------------------------------------
!   
      integer function Nr_mdle_sons(Type,Kref)
!
      implicit none
!
!  ...Arguments
      character(len=4), intent(in)  :: Type      
      integer,          intent(in)  :: Kref
!
      select case(Type)
!
!  ...brick
      case('mdlb')
        select case(Kref)
        case(111); Nr_mdle_sons = 8
        case(110,101,011); Nr_mdle_sons = 4
        case(100,010,001); Nr_mdle_sons = 2
        case default
          go to 999
        end select
!
!  ...prism
      case('mdlp')
        select case(Kref)
        case(11); Nr_mdle_sons = 8
        case(10); Nr_mdle_sons = 4
        case(01); Nr_mdle_sons = 2
        case default
          go to 999
        end select
!
!  ...tet
      case('mdln')
        select case(Kref)
        case(11,12,13); Nr_mdle_sons = 8
        case(24,32); Nr_mdle_sons = 4
        case default
          go to 999
        end select
!
!  ...pyramid
      case('mdld')
        select case(Kref)
        case(10); Nr_mdle_sons = 4
        case default
          go to 999
        end select
!
      case default
        go to 999
      end select
      return
!
!  ...error exit
  999 write(*,7100) Type,Kref
 7100 format('Nr_mdle_sons: Type,Kref = ',a4,3x,i3)
      stop 1
!
      end function Nr_mdle_sons
