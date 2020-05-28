!----------------------------------------------------------------------
!
!   routine name       - enforce_min_rule
!
!----------------------------------------------------------------------
!
!   latest revision    - Apr 18
!
!   purpose            - routine enforces the min rule for a FE mesh
!
!
!----------------------------------------------------------------------
!
   subroutine enforce_min_rule
!
   use data_structure3D
   use refinements
   use constrained_nodes
#include "syscom.blk"
!
!..work space for elem_nodes
   dimension nodesl(27),norientl(27)
!
!..order for element nodes implied by the order of the middle node
   dimension norder(19)
!
   character(len=4) :: etype
!
!----------------------------------------------------------------------
!
   iprint=0
!
!..reset visitation flags
   call reset_visit
!
!..loop through elements in the current mesh
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
      call get_connect_info(mdle, nodesl,norientl)
      call element_order(mdle,norientl, norder)
      if (iprint.eq.1) then
         write(*,7001) mdle,norder(1:nre+nrf+1)
 7001    format('enforce_min_rule: mdle = ',i5,' Norder = ',19i4)
      endif
!
!  ...loop through edges and faces
      do j=1,nre+nrf
         i=nrv+j
         nod = nodesl(i); nord = norder(j)
!
!     ...if active node
         if (NODES(nod)%act.eq.1) then
            call save_order(nod,nord)
!
!     ...constrained node
         else
!
!        ...identify the constraint case
            call decode2(NODES_CONSTR(i), nc,icase)
            select case(icase)
!
!        ...edge constrained by an edge
            case(11,12,37,38,47,48)
               nodp = NEDGC(nc)
               call save_order(nodp,nord)
!
!        ...horizontal edge constrained by a face
            case(26,28,33,36,63)
               nodp = NFACEC(nc)
               call save_order(nodp,nord*10+MAXP)
!
!        ...vertical edge constrained by a face
            case(25,27,43,46,53)
               nodp = NFACEC(nc)
               call save_order(nodp,MAXP*10+nord)
!
!        ...face node constrained by a face
            case(21,22,23,24,31,32,34,35,41,42,44,45,51,52,61,62)
               nodp = NFACEC(nc)
               call save_order(nodp,nord)
            end select
         endif
      enddo
   enddo
   if (iprint.eq.1) write(*,*) 'enforce_min_rule: DONE WITH Step 1'
!
!..loop through nodes
   do nod=1,NRNODS
!
!  ...active node
      if (NODES(nod)%act.eq.1) then
         nord = NODES(nod)%visit
         if (iprint.eq.1) write(*,*) 'enforce_min_rule: nod,nord = ',nod,nord
         if (nord.eq.0) cycle
!
!     ...refined node, check order of its sons
         if (NODES(nod)%ref_kind.ne.0) then
            select case(NODES(nod)%type)
!
!        ...edge node
            case('medg')
               do is=1,2
                  nods = NODES(nod)%sons(is)
                  nord = min(nord,NODES(nods)%visit)
               enddo
!
!        ...modify the node order
            if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
!        ...communicate the new order to the (inactive) sons
            do is=1,2
               nods = NODES(nod)%sons(is)
               NODES(nods)%order = nord
            enddo
!
!        ...triangular face node
            case('mdlt')
               call nr_face_sons('mdlt',NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  nods = NODES(nod)%sons(is)
                  select case(NODES(nods)%type)
                  case('mdlt')
                     nord = min(nord,NODES(nods)%visit)
                  case('mdlq')
                     call decode(NODES(nods)%visit, nordh,nordv)
                     if (nordh.ne.nordv) then
                        write(*,*) 'enforce_min_rule: INCONSISTENCY 3'
                        stop 1
                     endif
                     nord = min(nord,nordv)
                  end select
               enddo
!
!           ...modify the node order
               if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
!           ...communicate the new order to the (inactive) sons
               call nr_sons('mdlt',NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  nods = NODES(nod)%sons(is)
                  select case(NODES(nods)%type)
                  case('mdlt'); NODES(nods)%order = nord
                  case('mdlq'); NODES(nods)%order = nord*10+nord
                  end select
               enddo
!
!        ...rectangular face node
            case('mdlq')
               call decode(nord, nordh,nordv)
               call nr_face_sons('mdlq',NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  nods = NODES(nod)%sons(is)
                  call decode(NODES(nods)%visit, nordhs,nordvs)
                  nordh = min(nordh,nordhs); nordv = min(nordv,nordvs)
               enddo
!
!           ...modify the node order
               nord = nordh*10+nordv
               if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
!           ...communicate the new order to the (inactive) mid-face sons
               do is=1,nrsons
                  nods = NODES(nod)%sons(is)
                  NODES(nods)%order = nord
               enddo
!
!           ...communicate the new order to the (inactive) mid-edge sons
               select case(NODES(nod)%ref_kind)
               case(11)
                  do i=1,4
                     nods = NODES(nod)%sons(nrsons+i)
                     select case(i)
                     case(1,3); NODES(nods)%order = nordv
                     case(2,4); NODES(nods)%order = nordh
                     end select
                  enddo
               case(10)
                  nods = NODES(nod)%sons(nrsons+1)
                  NODES(nods)%order = nordv
               case(01)
                  nods = NODES(nod)%sons(nrsons+1)
                  NODES(nods)%order = nordh
               end select
            end select
!
!     ...unrefined active node
         else
!
!        ...modify the node order
            nord = NODES(nod)%visit
            if (iprint.eq.1) &
            write(*,*) 'enforce_min_rule: CALLING nodmod, nod,nord = ',nod,nord
            if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
!     ...if a refined node
         endif
!
!  ...if an active node
      endif
!
!..end of loop through nodes
   enddo
!
!..reset visitation flags
   call reset_visit
!
!
      end subroutine enforce_min_rule







   subroutine save_order(Nod,Nord)
   use data_structure3D
#include "syscom.blk"
!
   if (NODES(Nod)%visit.eq.0) then
      NODES(nod)%visit= Nord
   else
      select case(NODES(Nod)%type)
      case('medg','mdlt')
         NODES(nod)%visit = min(NODES(nod)%visit,Nord)
      case('mdlq')
         call decode(Nord, nordh1,nordv1)
         call decode(NODES(nod)%visit, nordh2,nordv2)
         nordh = min(nordh1,nordh2); nordv = min(nordv1,nordv2)
         NODES(nod)%visit = nordh*10+nordv
      end select
   endif
!
   end subroutine save_order
