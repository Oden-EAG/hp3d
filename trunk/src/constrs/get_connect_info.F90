!----------------------------------------------------------------------
!
!   routine name       - get_connect_info
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine returns element to nodes connectivities
!                        for an element, including info on constrained
!                        nodes
!
!   arguments :
!     in:
!             Mdle     - an element middle node, identified with
!                        the element
!     out:
!             Nodesl   - element nodes
!             Norientl - their orientations
!
!-----------------------------------------------------------------------
!
   subroutine get_connect_info(Mdle, Nodesl,Norientl)
!
      use element_data
      use data_structure3D
      use refinements
      use constrained_nodes, only: INFO_CONSTRAINTS,NODES_CONSTR, &
           NR_EDGES,NEDGC,NEDG_CONS,NR_FACES,NFACEC,NFACE_CONS,   &
           SON_NUM,FATH_TYPE,FATH_ORIENT,FATH_NODES,              &
           rotate_edge_nodes,rotate_trian_nodes,rotate_quadr_nodes
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Nodesl(27),Norientl(27)
!
      integer :: npar_refs(27),nson_refs(27),nort_refs(27)
!
!  ...miscellanea
      integer :: jef(4),jvf(4),kref_face(6)
      integer :: ntype,ftype,nfath
      integer :: i,ie,ifc,is,is1,j,je,jp,jv,jv1,jv2,loc
      integer :: iref,ireff,iref1,iref2,iref3,iref_fath
      integer :: nef,nff,nvf,nod,nodp,nort,nson,nrnodes
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!  ...initialize output
      Nodesl(1:27) = 0; Norientl(1:27) = 0
!
!  ...initialize variables like this or OpenMP will not work!
      ireff=0; is1=0;
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7000) Mdle
 7000   format('get_connect_info: Mdle = ',i6)
        call pause
      endif
#endif
!
      INFO_CONSTRAINTS=1
      call elem_nodes(Mdle, Nodesl,Norientl)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        call elem_show(Mdle, NODES(Mdle)%ntype,Nodesl,Norientl)
      endif
#endif
!
      INFO_CONSTRAINTS=0
!
!  ...store the element father
      nfath = NODES(Mdle)%father
!
!  ...initiate the data base for constraints
      NODES_CONSTR=0
      NR_EDGES=0; NEDGC=0; NEDG_CONS=0
      NR_FACES=0; NFACEC=0; NFACE_CONS=0
!
!  ...quit if an initial mesh element
      if (nfath.le.0) return
!
      nson = SON_NUM
      ftype = FATH_TYPE
      ntype = NODES(Mdle)%ntype
      iref_fath = NODES(nfath)%ref_kind
      call find_face_ref_flags(ftype,iref_fath, kref_face)
      call decode_ref(ftype,iref_fath, iref1,iref2,iref3)
!
!  ...number of vertices and edges for the father
      nvf = nvert(ftype); nef = nedge(ftype); nff = nface(ftype)
!
!  ...pre-compute info for loop
      call npar_ref_all(ftype,nson,iref1,iref2,iref3, npar_refs)
      call nson_ref_all(ftype,nson,iref1,iref2,iref3, nson_refs)
      call nort_ref_all(ftype,nson,iref1,iref2,iref3, nort_refs)
!
!  ...loop through nodes of the element son
      nrnodes = nvert(ntype)+nedge(ntype)+nface(ntype)+1
      do j=1,nrnodes
!
!  .....parent node
        jp   = npar_refs(j)
        is   = nson_refs(j)
        nort = nort_refs(j)
!
        if (is.ne.0) then
          nodp = FATH_NODES(jp)
!
!  .......parent edge node
          if (jp.le.nvf+nef) then
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,7021) j,jp,is
 7021         format('get_connect_info: 1: j,jp,is = ',3i5)
            endif
#endif
            call rotate_edge(FATH_ORIENT(jp),is,nort)
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,7022) j,jp,is
 7022         format('get_connect_info: 2: j,jp,is = ',3i5)
            endif
#endif
!            nod = NODES(nodp)%sons(is)
            nod = Son(nodp,is)
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,7023) nodp,nod,NODES(nod)%act
 7023         format('get_connect_info: nodp,nod,NODES(nod)%act = ', &
                     2i6,l2)
              call pause
            endif
#endif
!
!  .........if the node is constrained (inactive)
!!!            if (Is_inactive(nod)) then
            if (Is_inactive(nod) .and. NODES(nod)%ref_kind.eq.0) then
              call locate(nodp, NEDGC,NR_EDGES,loc)
              if (loc.eq.0) then
                NR_EDGES = NR_EDGES+1
                NEDGC(NR_EDGES) = nodp
                ie = jp-nvf
                call edge_to_vert(ftype,ie, jv1,jv2)
                NEDG_CONS(1,NR_EDGES) = FATH_NODES(jv1)
                NEDG_CONS(2,NR_EDGES) = FATH_NODES(jv2)
                call rotate_edge_nodes(FATH_ORIENT(jp),NR_EDGES)
                loc = NR_EDGES
              endif
              NODES_CONSTR(j) = loc*100+10+is
            endif
!
!  .......parent face node
          elseif (jp.le.nvf+nef+nff) then
!
!  .........local refinement flag for the face
            ifc = jp-nvf-nef
            iref = kref_face(ifc)
!
!  .........global refinement flag for the face
            ireff = NODES(nodp)%ref_kind
!
!  .........face edge and vertex nodes numbers
            call face_to_edge(ftype,ifc, jef(1),jef(2),jef(3),jef(4))
            call face_to_vert(ftype,ifc, jvf(1),jvf(2),jvf(3),jvf(4))
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,7100) j,jp,ifc,iref,ireff
 7100         format('get_connect_info: j,jp,ifc,iref,ireff = ',5i4)
            endif
#endif
!
            select case(TYPE_NOD(jp,ftype))
!
!  .........triangular parent face
            case(MDLT)
              call rotate_trian(iref,ireff,FATH_ORIENT(jp),is,nort)
!              nod = NODES(nodp)%sons(is)
              nod = Son(nodp,is)
#if HP3D_DEBUG
              if (iprint.eq.1) then
                write(*,7101) nodp,is,nod,NODES(nod)%act
 7101           format('get_connect_info: nodp,is,nod,NODES(nod)%act =', &
                                          3i5,l2)
              endif
#endif
              if (Is_inactive(nod).and.NODES(nod)%ref_kind.eq.0) then
!!!              if (Is_inactive(nod)) then
                call locate(nodp, NFACEC,NR_FACES,loc)
                if (loc.eq.0) then
                  NR_FACES = NR_FACES+1
                  NFACEC(NR_FACES) = nodp
                  do i=1,3
                    je = jef(i)
                    select case(FATH_ORIENT(nvf+je))
                    case (0)
                      NFACE_CONS(i,NR_FACES) = FATH_NODES(nvf+je)
                    case (1)
                      NFACE_CONS(i,NR_FACES) = -FATH_NODES(nvf+je)
                    end select
                    jv = jvf(i)
                    NFACE_CONS(4+i,NR_FACES) = FATH_NODES(jv)
                  enddo
                  call rotate_trian_nodes(FATH_ORIENT(jp),NR_FACES)
                  loc = NR_FACES
                endif
                select case(ireff)
                case(1)
                  NODES_CONSTR(j) = loc*100+70+is
                case(2,3,4)
                  select case(is)
                  case(1)
                    NODES_CONSTR(j) = loc*100+70+ireff-1
                  case(2)
                    NODES_CONSTR(j) = loc*100+80+ireff
                  case(3)
                    NODES_CONSTR(j) = loc*100+70+ireff+3
                  end select
                end select
              endif
!
!  .........quadrilateral parent face
            case(MDLQ)
              call rotate_quad(iref,ireff,FATH_ORIENT(jp), &
                               is,is1,nort)
!              nod = NODES(nodp)%sons(is)
              nod = Son(nodp,is)
              if (is1.ne.0) then
!                nod = NODES(nod)%sons(is1)
                nod = Son(nod,is1)
              endif
!
             if (Is_inactive(nod).and.NODES(nod)%ref_kind.eq.0) then
!!!             if (Is_inactive(nod)) then
                call locate(nodp, NFACEC,NR_FACES,loc)
                if (loc.eq.0) then
                  NR_FACES = NR_FACES+1
                  NFACEC(NR_FACES) = nodp
                  do i=1,4
                    je = jef(i)
                    select case(FATH_ORIENT(nvf+je))
                    case (0)
                      NFACE_CONS(i,NR_FACES) = FATH_NODES(nvf+je)
                    case (1)
                      NFACE_CONS(i,NR_FACES) = -FATH_NODES(nvf+je)
                    end select
                    jv = jvf(i)
                    NFACE_CONS(4+i,NR_FACES) = FATH_NODES(jv)
                  enddo
                  call rotate_quadr_nodes(FATH_ORIENT(jp),NR_FACES)
                  loc = NR_FACES
                endif
!
                select case(ireff)
                case(11)
                  NODES_CONSTR(j) = loc*100+20+is
                case(10)
                  if (iref.eq.11) then
                    call modify_face_info(nodp,is,is1, NODES_CONSTR(j))
                  else
                    NODES_CONSTR(j) = loc*100+50+is
                  endif
                case(01)
                  if (iref.eq.11) then
                    call modify_face_info(nodp,is,is1, NODES_CONSTR(j))
                  else
                    NODES_CONSTR(j) = loc*100+60+is
                  endif
                end select
              endif
            end select
          endif
        endif
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7011) Mdle
 7011   format('get_connect_info: CONSTRAINED NODES FOR Mdle = ',i6)
        write(*,7012) NODES_CONSTR(1:nrnodes)
 7012   format(27i5)
        write(*,7013)
 7013   format('                  CONSTRAINING EDGES = ')
        do i=1,NR_EDGES
          write(*,7014) NEDGC(i),NEDG_CONS(1:2,i)
 7014     format('EDGE = ',i6,' VERTICES = ',2i6)
        enddo
        write(*,7015)
 7015   format('                  CONSTRAINING FACES = ')
        do i=1,NR_FACES
          write(*,7016) NFACEC(i),NFACE_CONS(1:8,i)
 7016     format('FACE = ',i6,' EDGES = ',4i6,' VERTICES = ',4i6)
        enddo
        call pause
        call result
      endif
#endif
!
   end subroutine get_connect_info
!
!----------------------------------------------------------------------
!
!   routine name       - modify_face_info
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - a subsidiary of get_connect_info, routine
!                        modifies the data base for a constraining
!                        face that has been h2-refined
!
!   arguments :
!     in:
!                 Nodp - mid-face node
!                 Is   - son number for a son of the mid-face node
!                 Is1  - grandson number, and
!                        data base for constraining nodes in module
!                        'constrained_nodes'
!     out:
!                 Nodc - nickname for the constrained node, and
!                        modifications in the data base
!
!-----------------------------------------------------------------------
!
   subroutine modify_face_info(Nodp,Is,Is1, Nodc)
!
      use element_data
      use data_structure3D
      use refinements
      use constrained_nodes
!
      implicit none
!
      integer, intent(in)  :: Nodp,Is,Is1
      integer, intent(out) :: Nodc
!
      integer :: ifc,i,ie,ie1,iec,ise,iv,jv,loc,loc1
      integer :: medge,nson,nvoid,nvt
!
      integer, parameter :: ie_no(1:2,0:1) &
            = reshape( (/1,3, 4,2/) , (/2,2/) )
!
      integer, external :: imod
!!!      imod(j,mod) = j-(j-1)/mod*mod
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
!
      if (iprint.eq.1) then
        write(*,7001) Nodp,Is,Is1
 7001   format('modify_face_info: Nodp,Is,Is1 = ',i5,2i3)
      endif
#endif
!
      call locate(Nodp, NFACEC,NR_FACES,ifc)
      if (ifc.eq.0) then
        write(*,7002) Nodp
 7002   format('modify_face_info: INCONSISTENCY FOR Nodp = ',i5)
        stop 1
      endif
      select case(NODES(Nodp)%ref_kind)
!
!  ...the face has been h2-refined, vertically of horizontally
      case(10,01)
!
!  .....offset for the first edge
        call decode(NODES(Nodp)%ref_kind, nvoid,ie1)
!
!  .....move to the son
!        nson = NODES(Nodp)%sons(Is)
        nson = Son(Nodp,Is)
        select case(Is)
!
!  .....mid-face node son
        case(1,2)
!
!  .......set up a new constraining face
          call locate(nson,NFACEC,NR_FACES, loc)
          if (loc.eq.0) then
            NR_FACES = NR_FACES+1
!
!  .........store the mid-face node
            NFACEC(NR_FACES) = nson
            NFACE_CONS(1:8,NR_FACES)=0
!
!  .........store relevant edges only
            do ie=ie1+1,4,2
              medge = abs(NFACE_CONS(ie,ifc))
              if (NFACE_CONS(ie,ifc).gt.0) then
                ise = Is
              else
                ise = imod(Is+1,2)
              endif
              NFACE_CONS(ie,NR_FACES) = &
               sign(Son(medge,ise),NFACE_CONS(ie,ifc))
!              sign(NODES(medge)%sons(ise),NFACE_CONS(ie,ifc))
            enddo
            loc = NR_FACES
          endif
          Nodc = loc*100+30+ie1*10+(Is-1)*3+Is1
!
!  .....horizontal or vertical mid-edge node
        case(3)
!
!  .......set up a new constraining edge
          call locate(nson, NEDGC,NR_EDGES,iec)
          if (iec.eq.0) then
            NR_EDGES = NR_EDGES+1
            iec = NR_EDGES
            NEDGC(iec) = nson
            do iv=1,2
              ie = ie_no(iv,ie1)
              medge = abs(NFACE_CONS(ie,ifc))
!              nvt = NODES(medge)%sons(3)
              nvt = Son(medge,3)
!
!  ...........inactive vertex node, add the edge to the data base
              if (Is_inactive(nvt)) then
                call locate(medge, NEDGC,NR_EDGES,loc1)
                if (loc1.eq.0) then
                  NR_EDGES = NR_EDGES+1
                  NEDGC(NR_EDGES) = medge
                  do i=1,2
                    jv = QUADR_EDGE_TO_VERT(i,ie)
                    NEDG_CONS(i,NR_EDGES) = NFACE_CONS(4+jv,ifc)
                  enddo
                  if (NFACE_CONS(ie,ifc).lt.0) &
                    call rotate_edge_nodes(1,NR_EDGES)
                  loc1 = NR_EDGES
                endif
                nvt = -(loc1*100+13)
              endif
              NEDG_CONS(iv,iec) = nvt
            enddo
          endif
          Nodc = iec*100+ie1*10+36+Is1
        end select
!
!  ...any other option illegal
      case default
        write(*,7003) Nodp
 7003   format('modify_face_info: Nodp = ',i5)
        stop 1
      end select
!
!
   end subroutine modify_face_info
