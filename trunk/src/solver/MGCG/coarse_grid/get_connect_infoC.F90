!----------------------------------------------------------------------
!
!   routine name       - get_connect_infoC
!
!----------------------------------------------------------------------
!
!   latest revision    - Jan 18
!
!   purpose            - routine returns element to nodes connectivities
!                        for a coarse grid element, including info on
!                        constrained nodes
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
   subroutine get_connect_infoC(Igrid,Mdle, Nodesl,Norientl)
!
   use element_data
   use data_structure3D
   use refinements,       only: decode_ref, npar_ref, nson_ref, &
                                nort_ref, rotate_edge, rotate_trian,  &
                                rotate_quad, elem_show
   use constrained_nodes, only: INFO_CONSTRAINTS,NODES_CONSTR, &
                                NR_EDGES,NEDGC,NEDG_CONS,NR_FACES,NFACEC,  &
                                NFACE_CONS, SON_NUM,FATH_TYPE,FATH_ORIENT, &
                                FATH_NODES, rotate_edge_nodes,rotate_trian_nodes,  &
                                rotate_quadr_nodes
   use mg_data_structure
#include "syscom.blk"
!
   dimension Nodesl(27),Norientl(27)
!
!..face edge and vertex numbers etc
   integer, intent(in) :: Igrid
   character(len=4) :: type,ftype,stype
   integer :: if, ireff, is1, jef(4), jvf(4), kref_face(6)
!
!..initialize variables like this or OpenMP will not work!
   ireff=0; is1=0;
!
   select case(Mdle)
   case(135)
      iprint=0
   case default
      iprint=0
   end select
   if (iprint.eq.1) then
      write(*,7000) Mdle
 7000 format('get_connect_info: Mdle = ',i6)
      call pause
   endif
!
!
   INFO_CONSTRAINTS=1
   call elem_nodes(Mdle, Nodesl,Norientl)
   if (iprint.eq.1) then
      call elem_show(Mdle, NODES(Mdle)%type,Nodesl,Norientl)
   endif
   INFO_CONSTRAINTS=0
!
!..store the element father
   nfath = NODES(Mdle)%father
!
!..initiate the data base for constraints
   NODES_CONSTR=0
   NR_EDGES=0; NEDGC=0; NEDG_CONS=0
   NR_FACES=0; NFACEC=0; NFACE_CONS=0
!
!..quit if an initial mesh element
   if (nfath.le.0) return
!
   nson = SON_NUM
   ftype = FATH_TYPE
   stype = NODES(Mdle)%type
   iref_fath = NODES(nfath)%ref_kind
   call find_face_ref_flags(ftype,iref_fath, kref_face)
   call decode_ref(ftype, iref_fath, iref1, iref2, iref3)
!
!..number of vertices and edges for the father
   nvf = nvert(ftype); nef = nedge(ftype); nff = nface(ftype)
!
!..loop through nodes of the element son
   do j=1,nvert(stype)+nedge(stype)+nface(stype)+1
!
!  ...parent node
      jp   = npar_ref(ftype, j, nson, iref1,iref2,iref3)
      is   = nson_ref(ftype, j, nson, iref1,iref2,iref3)
      nort = nort_ref(ftype, j, nson, iref1,iref2,iref3)
      if (is.ne.0) then
         nodp = FATH_NODES(jp)
!
!     ...parent edge node
         if (jp.le.nvf+nef) then
            if (iprint.eq.1) then
               write(*,7021) j,jp,is
 7021          format('get_connect_info: 1: j,jp,is = ',3i5)
            endif
            call rotate_edge(FATH_ORIENT(jp),is,nort)
            if (iprint.eq.1) then
               write(*,7022) j,jp,is
 7022          format('get_connect_info: 2: j,jp,is = ',3i5)
            endif
!            nod = NODES(nodp)%sons(is)
            nod = Son(nodp,is)
            if (iprint.eq.1) then
               write(*,7023) nodp,nod,NODES(nod)%act
 7023          format('get_connect_info: nodp,nod,NODES(nod)%act = ' &
                      2i6,i2)
               call pause
            endif
!
!        ...if the coarse grid node was constrained (inactive), its
!           parent node has been marked as an active  node
            if (NODES_MG(nodp)%master(Igrid).eq.1) then
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
!     ...parent face node
         elseif (jp.le.nvf+nef+nff) then
!
!        ...local refinement flag for the face
            if = jp-nvf-nef
            iref = kref_face(if)
!
!        ...global refinement flag for the face
            ireff = NODES(nodp)%ref_kind
!
!        ...face edge and vertex nodes numbers
            call face_to_edge(ftype,if, jef(1),jef(2),jef(3),jef(4))
            call face_to_vert(ftype,if, jvf(1),jvf(2),jvf(3),jvf(4))
            if (iprint.eq.1) then
               write(*,7100) j,jp,if,iref,ireff
 7100          format('get_connect_info: j,jp,if,iref,ireff = ',5i4)
            endif
!
            select case(type_nod(ftype,jp))
!
!        ...triangular parent face
            case('mdlt')
               call rotate_trian(iref,ireff,FATH_ORIENT(jp),is,nort)
!               nod = NODES(nodp)%sons(is)
               nod = Son(nodp,is)
               if (iprint.eq.1) then
                  write(*,7101) nodp,is,nod,NODES(nod)%act
 7101             format('get_connect_info: nodp,is,nod,NODES(nod)%act =',4i5)
               endif
!
!
!          ...if the coarse grid node was constrained (inactive), its
!             parent node has been marked as an active  node
               if (NODES_MG(nodp)%master(Igrid).eq.1) then
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
!        ...quadrilateral parent face
            case('mdlq')
               call rotate_quad(iref,ireff,FATH_ORIENT(jp), is,is1,nort)
!               nod = NODES(nodp)%sons(is)
               nod = Son(nodp,is)
               if (is1.ne.0) then
!                  nod = NODES(nod)%sons(is1)
                  nod = Son(nod,is1)
               endif
!
!           ...if the coarse grid node was constrained (inactive), its
!              parent node has been marked as an active  node
               if (NODES_MG(nodp)%master(Igrid).eq.1) then
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
   if (iprint.eq.1) then
      write(*,7011) Mdle
 7011 format('get_connect_infoC: CONSTRAINED NODES FOR Mdle = ',i6)
      type = NODES(Mdle)%type
      nrn = nvert(type)+nedge(type)+nface(type)+1
      write(*,7012) NODES_CONSTR(1:nrn)
 7012 format(27i5)
      write(*,7013)
 7013 format('                  CONSTRAINING EDGES = ')
      do i=1,NR_EDGES
         write(*,7014) NEDGC(i),NEDG_CONS(1:2,i)
 7014    format('EDGE = ',i6,' VERTICES = ',2i6)
      enddo
      write(*,7015)
 7015 format('                  CONSTRAINING FACES = ')
      do i=1,NR_FACES
         write(*,7016) NFACEC(i),NFACE_CONS(1:8,i)
 7016    format('FACE = ',i6,' EDGES = ',4i6,' VERTICES = ',4i6)
      enddo
         call pause
         call result
      endif
!
      end subroutine get_connect_infoC
