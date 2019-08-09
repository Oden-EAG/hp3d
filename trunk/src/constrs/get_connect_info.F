c----------------------------------------------------------------------
c
c   routine name       - get_connect_info
c
c----------------------------------------------------------------------
c
c   latest revision    - Jul 10
c
c   purpose            - routine returns element to nodes connectivities
c                        for an element, including info on constrained
c                        nodes
c
c   arguments :
c     in:
c             Mdle     - an element middle node, identified with
c                        the element
c     out:
c             Nodesl   - element nodes
c             Norientl - their orientations
c
c-----------------------------------------------------------------------
c
      subroutine get_connect_info(Mdle, Nodesl,Norientl)
c
      use element_data
      use data_structure3D
      use refinements, only : decode_ref, npar_ref, nson_ref,
     .    nort_ref, rotate_edge, rotate_trian, rotate_quad, elem_show
      use constrained_nodes, only: INFO_CONSTRAINTS,NODES_CONSTR,
     .     NR_EDGES,NEDGC,NEDG_CONS,NR_FACES,NFACEC,NFACE_CONS,
     .     SON_NUM,FATH_TYPE,FATH_ORIENT,FATH_NODES,
     .     rotate_edge_nodes,rotate_trian_nodes,rotate_quadr_nodes
#include "syscom.blk"
c
      dimension Nodesl(27),Norientl(27)
c
c  ...face edge and vertex numbers etc
      character(len=4) :: type,ftype,stype
      integer :: if, ireff, is1, jef(4), jvf(4), kref_face(6)
c
c  ...initialize variables like this or OpenMP will not work!
      ireff=0; is1=0;
c
      select case(Mdle)
      case(135)
        iprint=0
      case default
        iprint=0
      end select
      if (iprint.eq.1) then
        write(*,7000) Mdle
 7000   format('get_connect_info: Mdle = ',i6)
        call pause
      endif
c
c
      INFO_CONSTRAINTS=1
      call elem_nodes(Mdle, Nodesl,Norientl)
      if (iprint.eq.1) then
        call elem_show(Mdle, NODES(Mdle)%type,Nodesl,Norientl)
      endif
      INFO_CONSTRAINTS=0
c
c  ...store the element father
      nfath = NODES(Mdle)%father
c
c  ...initiate the data base for constraints
      NODES_CONSTR=0
      NR_EDGES=0; NEDGC=0; NEDG_CONS=0
      NR_FACES=0; NFACEC=0; NFACE_CONS=0
c
c  ...quit if an initial mesh element
      if (nfath.le.0) return
c
      nson = SON_NUM
      ftype = FATH_TYPE
      stype = NODES(Mdle)%type
      iref_fath = NODES(nfath)%ref_kind
      call find_face_ref_flags(ftype,iref_fath, kref_face)
      call decode_ref(ftype, iref_fath, iref1, iref2, iref3)
c
c  ...number of vertices and edges for the father
      nvf = nvert(ftype); nef = nedge(ftype); nff = nface(ftype)
c
c  ...loop through nodes of the element son
      do j=1,nvert(stype)+nedge(stype)+nface(stype)+1
c
c  .....parent node
        jp   = npar_ref(ftype, j, nson, iref1,iref2,iref3)
        is   = nson_ref(ftype, j, nson, iref1,iref2,iref3)
        nort = nort_ref(ftype, j, nson, iref1,iref2,iref3)
        if (is.ne.0) then
          nodp = FATH_NODES(jp)
c
c  .......parent edge node
          if (jp.le.nvf+nef) then
            if (iprint.eq.1) then
              write(*,7021) j,jp,is
 7021         format('get_connect_info: 1: j,jp,is = ',3i5)
            endif
            call rotate_edge(FATH_ORIENT(jp),is,nort)
            if (iprint.eq.1) then
              write(*,7022) j,jp,is
 7022         format('get_connect_info: 2: j,jp,is = ',3i5)
            endif
c            nod = NODES(nodp)%sons(is)
            nod = Son(nodp,is)
            if (iprint.eq.1) then
              write(*,7023) nodp,nod,NODES(nod)%act
 7023         format('get_connect_info: nodp,nod,NODES(nod)%act = ',
     .               2i6,i2)
              call pause
            endif
c
c  .........if the node is constrained (inactive)
ccc            if (NODES(nod)%act.eq.0) then
            if (NODES(nod)%act.eq.0 .and. NODES(nod)%ref_kind.eq.0) then
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
c
c  .......parent face node
          elseif (jp.le.nvf+nef+nff) then
c
c  .........local refinement flag for the face
            if = jp-nvf-nef
            iref = kref_face(if)
c
c  .........global refinement flag for the face
            ireff = NODES(nodp)%ref_kind
c
c  .........face edge and vertex nodes numbers
            call face_to_edge(ftype,if, jef(1),jef(2),jef(3),jef(4))
            call face_to_vert(ftype,if, jvf(1),jvf(2),jvf(3),jvf(4))
            if (iprint.eq.1) then
              write(*,7100) j,jp,if,iref,ireff
 7100         format('get_connect_info: j,jp,if,iref,ireff = ',5i4)
            endif
c
            select case(type_nod(ftype,jp))
c
c  .........triangular parent face
            case('mdlt')
              call rotate_trian(iref,ireff,FATH_ORIENT(jp),is,nort)
c              nod = NODES(nodp)%sons(is)
              nod = Son(nodp,is)
              if (iprint.eq.1) then
                write(*,7101) nodp,is,nod,NODES(nod)%act
 7101           format('get_connect_info: nodp,is,nod,NODES(nod)%act =',
     .                                    4i5)
              endif
              if (NODES(nod)%act.eq.0.and.NODES(nod)%ref_kind.eq.0) then
ccc              if (NODES(nod)%act.eq.0) then
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
c
c  .........quadrilateral parent face
            case('mdlq')
              call rotate_quad(iref,ireff,FATH_ORIENT(jp),
     .                         is,is1,nort)
c              nod = NODES(nodp)%sons(is)
              nod = Son(nodp,is)
              if (is1.ne.0) then
c                nod = NODES(nod)%sons(is1)
                nod = Son(nod,is1)
              endif
c
             if (NODES(nod)%act.eq.0.and.NODES(nod)%ref_kind.eq.0) then
ccc             if (NODES(nod)%act.eq.0) then
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
c
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
c
      if (iprint.eq.1) then
        write(*,7011) Mdle
 7011   format('get_connect_info: CONSTRAINED NODES FOR Mdle = ',i6)
        type = NODES(Mdle)%type
        nrn = nvert(type)+nedge(type)+nface(type)+1
        write(*,7012) NODES_CONSTR(1:nrn)
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
c
      end subroutine
c
c----------------------------------------------------------------------
c
c   routine name       - modify_face_info
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - a subsidiary of get_connect_info, routine
c                        modifies the data base for a constraining
c                        face that has been h2-refined
c
c   arguments :
c     in:
c                 Nodp - mid-face node
c                 Is   - son number for a son of the mid-face node
c                 Is1  - grandson number, and
c                        data base for constraining nodes in module
c                        'constrained_nodes'
c     out:
c                 Nodc - nickname for the constrained node, and
c                        modifications in the data base
c
c-----------------------------------------------------------------------
c
      subroutine modify_face_info(Nodp,Is,Is1, Nodc)
c
      use element_data
      use data_structure3D
      use refinements
      use constrained_nodes
#include "syscom.blk"
c
      integer, parameter, dimension(1:2,0:1) :: ie_no
     .  = reshape( (/1,3, 4,2/) , (/2,2/) )
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
c
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) Nodp,Is,Is1
 7001   format('modify_face_info: Nodp,Is,Is1 = ',i5,2i3)
      endif
c
      call locate(Nodp, NFACEC,NR_FACES,ifc)
      if (ifc.eq.0) then
        write(*,7002) Nodp
 7002   format('modify_face_info: INCONSISTENCY FOR Nodp = ',i5)
        stop 1
      endif
      select case(NODES(Nodp)%ref_kind)
c
c  ...the face has been h2-refined, vertically of horizontally
      case(10,01)
c
c  .....offset for the first edge
        call decode(NODES(Nodp)%ref_kind, nvoid,ie1)
c
c  .....move to the son
c        nson = NODES(Nodp)%sons(Is)
        nson = Son(Nodp,Is)
        select case(Is)
c
c  .....mid-face node son
        case(1,2)
c
c  .......set up a new constraining face
          call locate(nson,NFACEC,NR_FACES, loc)
          if (loc.eq.0) then
            NR_FACES = NR_FACES+1
c
c  .........store the mid-face node
            NFACEC(NR_FACES) = nson
            NFACE_CONS(1:8,NR_FACES)=0
c
c  .........store relevant edges only
            do ie=ie1+1,4,2
              medg = iabs(NFACE_CONS(ie,ifc))
              if (NFACE_CONS(ie,ifc).gt.0) then
                ise = Is
              else
                ise = imod(Is+1,2)
              endif
              NFACE_CONS(ie,NR_FACES) =
     .        sign(Son(medg,ise),NFACE_CONS(ie,ifc))
c     .        sign(NODES(medg)%sons(ise),NFACE_CONS(ie,ifc))
            enddo
            loc = NR_FACES
          endif
          Nodc = loc*100+30+ie1*10+(Is-1)*3+Is1
c
c  .....horizontal or vertical mid-edge node
        case(3)
c
c  .......set up a new constraining edge
          call locate(nson, NEDGC,NR_EDGES,iec)
          if (iec.eq.0) then
            NR_EDGES = NR_EDGES+1
            iec = NR_EDGES
            NEDGC(iec) = nson
            do iv=1,2
              ie = ie_no(iv,ie1)
              medg = iabs(NFACE_CONS(ie,ifc))
c              nvt = NODES(medg)%sons(3)
              nvt = Son(medg,3)
c
c  ...........inactive vertex node, add the edge to the data base
              if (NODES(nvt)%act.eq.0) then
                call locate(medg, NEDGC,NR_EDGES,loc1)
                if (loc1.eq.0) then
                  NR_EDGES = NR_EDGES+1
                  NEDGC(NR_EDGES) = medg
                  do i=1,2
                    jv = QUADR_EDGE_TO_VERT(i,ie)
                    NEDG_CONS(i,NR_EDGES) = NFACE_CONS(4+jv,ifc)
                  enddo
                  if (NFACE_CONS(ie,ifc).lt.0)
     .              call rotate_edge_nodes(1,NR_EDGES)
                  loc1 = NR_EDGES
                endif
                nvt = -(loc1*100+13)
              endif
              NEDG_CONS(iv,iec) = nvt
            enddo
          endif
          Nodc = iec*100+ie1*10+36+Is1
        end select
c
c  ...any other option illegal
      case default
        write(*,7003) Nodp
 7003   format('modify_face_info: Nodp = ',i5)
        stop 1
      end select
c
c
      end subroutine


