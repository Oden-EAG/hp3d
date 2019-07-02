c----------------------------------------------------------------------
c
c   routine name       - elem_nodes
c
c----------------------------------------------------------------------
c
c   latest revision    - Jun 10
c
c   purpose            - routine returns nodes for an element (active
c                        or not) and their orientations
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
      subroutine elem_nodes(Mdle, Nodesl,Norientl)
c
      use element_data
      use data_structure3D
      use refinements
      use constrained_nodes
#include "syscom.blk"
c
      dimension Nodesl(27),Norientl(27)
c
c  ...history information: father, its type, refinement kind, son number
      character(len=4) :: stype,ftype,father_type,son_type
      dimension nfather(maxgen),father_type(maxgen),son_type(maxgen),
     .          nfather_ref_kind(maxgen),noson(maxgen)
c
c  ...local copies of Nodesl,Norientl
      dimension nodesl_loc(27),norientl_loc(27)
      integer :: ireff=0, is1=0
c
c  ...local refinement flags for an element faces
      dimension kref_face(6)
c
      select case(Mdle)
      case(155)
        iprint=0
      case default
        iprint=0
      end select
c
      if (iprint.eq.1) then
        write(*,*) '------------------------------------------------'
        write(*,7001) Mdle
 7001   format('elem_nodes: DEBUGGING FOR Mdle = ',i6)
      endif
c
      Nodesl=0; Norientl=0
c
c  ...Step 1: Go up the tree to the initial mesh ancestor
c     initiate father, son and generation number
      nfath = NODES(Mdle)%father; nson = Mdle; igen=0
      do while(nfath.gt.0)
        igen=igen+1
        nfather(igen) = nfath
        father_type(igen) = NODES(nfath)%type
        son_type(igen) = NODES(nson)%type
        nfather_ref_kind(igen) = NODES(nfath)%ref_kind
        call nr_mdle_sons(NODES(nfath)%type,NODES(nfath)%ref_kind,
     .                    nrsons)
        call locate(nson,NODES(nfath)%sons,nrsons, noson(igen))
        nson = nfath
        nfath = NODES(nson)%father
      enddo
      nrgen = igen
c
c  ...the initial mesh ancestor is "nson"
      nel = nson
      if (iprint.eq.1) then
        do igen=1,nrgen
          write(*,7011) igen,nfather(igen),father_type(igen),
     .    son_type(igen),nfather_ref_kind(igen),noson(igen)
 7011     format('igen = ',i2,' father = ',i6,' type = ',a5,
     .           ' son type = ',a5,
     .           ' ref_kind = ',i3,' son number = ',i1)
        enddo
      endif
c
c  ...collect the info for the initial mesh ancestor from ELEMS
      select case(ELEMS(nel)%type)
      case('tetr')
        Nodesl(1:15) = ELEMS(nel)%nodes(1:15)
        call decodg(ELEMS(nel)%edge_orient,2,6, Norientl(5:10))
        call decodg(ELEMS(nel)%face_orient,8,4, Norientl(11:14))
      case('pris')
        Nodesl(1:21) = ELEMS(nel)%nodes(1:21)
        call decodg(ELEMS(nel)%edge_orient,2,9, Norientl(7:15))
        call decodg(ELEMS(nel)%face_orient,8,5, Norientl(16:20))
      case('pyra')
        Nodesl(1:19) = ELEMS(nel)%nodes(1:19)
        call decodg(ELEMS(nel)%edge_orient,2,8, Norientl(6:13))
        call decodg(ELEMS(nel)%face_orient,8,5, Norientl(14:18))
      case('bric')
        Nodesl(1:27) = ELEMS(nel)%nodes(1:27)
        call decodg(ELEMS(nel)%edge_orient,2,12, Norientl(9:20))
        call decodg(ELEMS(nel)%face_orient,8,6, Norientl(21:26))
      end select
c
      if (iprint.eq.1) then
        write(*,7021) nel
 7021   format('elem_nodes: NODES AND ORIENTATIONS FOR nel = ',i6)
        select case(ELEMS(nel)%type)
        case('tetr')
          write(*,7103) Nodesl(1:15)
          write(*,7103) Norientl(1:15)
        case('pris')
          write(*,7104) Nodesl(1:21)
          write(*,7104) Norientl(1:21)
        case('pyra')
          write(*,7105) Nodesl(1:19)
          write(*,7105) Norientl(1:19)
        case('bric')
          write(*,7106) Nodesl(1:27)
          write(*,7106) Norientl(1:27)
        end select
      endif
c
c  ...Step 2: Go down the tree reconstructing connectivities
      if (nrgen.eq.0) return
      do igen=nrgen,1,-1
        nson = noson(igen)
        stype = son_type(igen); ftype = father_type(igen)
        iref_fath = nfather_ref_kind(igen)
        call find_face_ref_flags(ftype,iref_fath, kref_face)
        call decode_ref(ftype, iref_fath, iref1, iref2, iref3)
c
c  .....loop through nodes of the son
        do j=1,nvert(stype)+nedge(stype)+nface(stype)+1
c
c  .......parent node
          jp   = npar_ref(ftype, j, nson, iref1,iref2,iref3)
          is   = nson_ref(ftype, j, nson, iref1,iref2,iref3)
          nort = nort_ref(ftype, j, nson, iref1,iref2,iref3)

          if (iprint.eq.1) then
            write(*,7027) stype,j,nson,ftype,jp,is,nort
 7027       format('elem_nodes: stype,j,nson,ftype,jp,is,nort = ',
     .              a5,2i3,2x,a5,3i3)
            call pause
          endif
c
c  .......a node shared with the father
          if (is.eq.0) then
c  .........Kyungjoo temporary
            if (jp.eq.0) then
              nodesl_loc(j) = Nodesl(j)
              norientl_loc(j) = Norientl(j)
            else
              nodesl_loc(j) = Nodesl(jp)
              norientl_loc(j) = Norientl(jp)
            endif
c
c  .......a node resulted from breaking
          else
c
c  .........parent node
            nodp = Nodesl(jp)
c
            if (jp.le.nvert(ftype)) then
              write(*,7999) 1; stop 1
c
c  .........parent edge node
            elseif (jp.le.nvert(ftype)+nedge(ftype)) then
              call rotate_edge(Norientl(jp),is,nort)
              nodesl_loc(j) = NODES(nodp)%sons(is)
c
c  .........parent face node
            elseif (jp.le.nvert(ftype)+nedge(ftype)+nface(ftype)) then
c
c  ...........local refinement flag for the face
              iref = kref_face(jp-nvert(ftype)-nedge(ftype))
c
c  ...........global refinement flag for the face
              ireff = NODES(nodp)%ref_kind
c
              select case(type_nod(ftype, jp))
c
c  ...........triangular parent face
              case('mdlt')
                call rotate_trian(iref,ireff,Norientl(jp),is,nort)
                nodesl_loc(j) = NODES(nodp)%sons(is)
c
c  ...........quadrilateral parent face
              case('mdlq')
                call rotate_quad(iref,ireff,Norientl(jp),
     .                           is,is1,nort)
                if (iprint.eq.1) then
                  write(*,7039) iref,ireff,Norientl(jp),is,is1,nort
 7039             format('elem_nodes: AFTER ROTATION',
     .                   ' iref,ireff,Norientl(jp),is,is1,nort = ',
     .                    2i3,2x,i2,3x,2i2,3x,i2)
                  write(*,7041) is,nodp,NODES(nodp)%sons(is)
 7041             format('TAKING ',i2,'th SON OF ',i6,' = ',i6)
                endif
                nodesl_loc(j) = NODES(nodp)%sons(is)
                if (is1.ne.0) then
                  nodpp = nodesl_loc(j)
                  nodesl_loc(j) = NODES(nodpp)%sons(is1)
                  if (iprint.eq.1) then
                    write(*,7041) is1,nodpp,NODES(nodpp)%sons(is1)
                  endif
                endif
              end select
c
c  .........parent middle node
            else
              nodesl_loc(j) = NODES(nodp)%sons(is)
            endif
            norientl_loc(j) = nort
c
c  .......if the node has resulted from breaking a node
          endif
c
c  .....end of loop through nodes of the son
        enddo
c
        if (igen.eq.1) then
          if (INFO_CONSTRAINTS.eq.1) then
c
c  .........save the data for collecting info on constraints
            FATH_NODES = Nodesl
            FATH_ORIENT = Norientl
            FATH_TYPE = father_type(1)
            SON_NUM = nson
          endif
        endif
        Nodesl = nodesl_loc; Norientl = norientl_loc
        if (iprint.eq.1) then
          write(*,7031) igen,NODES(nfather(igen))%sons(nson)
 7031     format('elem_nodes: igen = ',i2,
     .           ' NODES AND ORIENT FOR mdle = ',i6)
          select case(son_type(igen))
          case('mdln')
            write(*,7103) Nodesl(1:15)
            write(*,7103) Norientl(1:15)
          case('mdlp')
            write(*,7104) Nodesl(1:21)
            write(*,7104) Norientl(1:21)
          case('mdld')
            write(*,7105) Nodesl(1:19)
            write(*,7105) Norientl(1:19)
          case('mdlb')
            write(*,7106) Nodesl(1:27)
            write(*,7106) Norientl(1:27)
          end select
        endif
      enddo
      if (iprint.eq.1) call pause
c
 7999 format('elem_nodes: INCONSISTENCY ',i2)
c
c  ...formats for printing element nodes and orientations
 7103 format(4i6,2x,6i6,2x,4i6,2x,i6)
 7104 format(6i6,2x,9i6,2x,2i6,2x,3i6,2x,i6)
 7105 format(5i6,2x,8i6,2x,i6,2x,4i6,2x,i6)
 7106 format(8i6,2x,12i6,2x,6i6,2x,i6)
c
      if ((iprint.eq.2).and.(INFO_CONSTRAINTS.eq.1)) then
        write(*,7200) Mdle
 7200   format('elem_nodes: FATHER INFO FOR Mdle = ',i6)
        nrs = nvert(FATH_TYPE)+nedge(FATH_TYPE)+nface(FATH_TYPE)+1
        write(*,7201) FATH_NODES(1:nrs)
 7201   format('FATHER NODES = ',27i6)
        write(*,7202) FATH_ORIENT(1:nrs)
 7202   format('NODES ORIENT = ',27i6)
        write(*,7203) FATH_TYPE,SON_NUM
 7203   format('FATH_TYPE = ',a5,' SON_NUM = ',i1)
        call pause
      endif
c
      end subroutine
