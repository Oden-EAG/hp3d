!----------------------------------------------------------------------
!
!   routine name       - enforce_min_rule
!
!----------------------------------------------------------------------
!
!   latest revision    - Aug 17
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
#include "syscom.blk"
!
!  ...work space for elem_nodes
      dimension nodesl(27),norientl(27)
!
!  ...order for element nodes implied by the order of the middle node
      dimension norder(19)
!
      character(len=4) :: etype
!
!----------------------------------------------------------------------
!
      iprint=0
!
!  ...reset visitation flags
      call reset_visit
!
!  ...loop through elements in the current mesh
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        etype = NODES(Mdle)%type
        nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
        call elem_nodes(mdle, nodesl,norientl)
        call element_order(mdle,norientl, norder)
        if (iprint.eq.1) then
          write(*,7001) mdle,norder(1:nre+nrf+1)
 7001     format('enforce_min_rule: mdle = ',i5,' Norder = ',19i4)
        endif
!     ..minimum rule for faces and edges wrt elements mdle node order       
        do j=1,nre+nrf
          nod = nodesl(nrv+j)
!
!  .......save the order implied
          if (NODES(nod)%visit.eq.0) then
            NODES(nod)%visit= norder(j)
          else
            select case(NODES(nod)%type)
            case('medg','mdlt')
              NODES(nod)%visit = min(NODES(nod)%visit,norder(j))
            case('mdlq')
              call decode(norder(j), nordh1,nordv1)
              call decode(NODES(nod)%visit, nordh2,nordv2)
              nordh = min(nordh1,nordh2); nordv = min(nordv1,nordv2)
              NODES(nod)%visit = nordh*10+nordv
            end select
          endif
        enddo 
      enddo
      if (iprint.eq.1) write(*,*) 'enforce_min_rule: DONE WITH Step 1'
!

!  ...loop through nodes
      do nod=1,NRNODS
        if (NODES(nod)%type .ne. 'mdlq' .and. NODES(nod)%type .ne. 'mdlt') cycle

        if (nod .eq. 1504 ) then 
          write(*,*) 'nod, NODES(nod)%visit=',nod, NODES(nod)%visit
        endif  
!
!  .....active node
        if (NODES(nod)%act.eq.1) then
          nord = NODES(nod)%visit
          if (iprint.eq.1) write(*,*) 'enforce_min_rule: nod,nord = ',nod,nord
          if (nord.eq.0) cycle
!
!  .......refined node, check order of its sons
          if (NODES(nod)%ref_kind.ne.0) then
            select case(NODES(nod)%type)
!
!  .........triangular face node
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
                if (NODES(nods)%ref_kind.ne.0) then
                  write(*,*) 'enforce_min_rule: INCONSISTENCY 4'
                  stop 1
                endif
              enddo
!
!  ...........modify the node order
              if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
!  ...........communicate the new order to the (inactive) sons
              call nr_sons('mdlt',NODES(nod)%ref_kind, nrsons)
              do is=1,nrsons
                nods = NODES(nod)%sons(is)
                if (NODES(nods)%act.ne.0) then
                  write(*,*) 'enforce_min_rule: INCONSISTENCY 5'
                  stop 1
                endif
                select case(NODES(nods)%type)
                case('mdlt'); NODES(nods)%order = nord
                case('mdlq'); NODES(nods)%order = nord*10+nord
                end select
              enddo
!
!  .........rectangular face node
            case('mdlq')
              call decode(nord, nordh,nordv)
              call nr_face_sons('mdlq',NODES(nod)%ref_kind, nrsons)
              do is=1,nrsons
                nods = NODES(nod)%sons(is)
                call decode(NODES(nods)%visit, nordhs,nordvs)
                nordh = min(nordh,nordhs); nordv = min(nordv,nordvs)
                if (NODES(nods)%ref_kind.ne.0) then
                  write(*,*) 'enforce_min_rule: INCONSISTENCY 6'
                  stop 1
                endif
              enddo
!
!  ...........modify the node order
              nord = nordh*10+nordv
              if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
!  ...........communicate the new order to the (inactive) mid-face sons
              do is=1,nrsons
                nods = NODES(nod)%sons(is)
                if (NODES(nods)%act.ne.0) then
                  write(*,*) 'enforce_min_rule: INCONSISTENCY 7'
                  stop 1
                endif
                NODES(nods)%order = nord
              enddo
!
!  ...........communicate the new order to the (inactive) mid-edge sons
              select case(NODES(nod)%ref_kind)
              case(11)
                do i=1,4
                  nods = NODES(nod)%sons(nrsons+i)
                  select case(i)
                  case(1,3); NODES(nods)%order = nordv
                  case(2,4); NODES(nods)%order = nordh
                  end select
                  if (NODES(nods)%ref_kind.ne.0) then
                    write(*,*) 'enforce_min_rule: INCONSISTENCY 8'
                    stop 1
                  endif
                enddo
              case(10)
                nods = NODES(nod)%sons(nrsons+1)
                NODES(nods)%order = nordv
                if (NODES(nods)%ref_kind.ne.0) then
                  write(*,*) 'enforce_min_rule: INCONSISTENCY 9'
                  stop 1
                endif
              case(01)
                nods = NODES(nod)%sons(nrsons+1)
                NODES(nods)%order = nordh
                if (NODES(nods)%ref_kind.ne.0) then
                  write(*,*) 'enforce_min_rule: INCONSISTENCY 10'
                  stop 1
                endif
              end select
            end select
!
!  .......unrefined active node
          else
!
!  .........modify the node order
            nord = NODES(nod)%visit
            if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
!  .......if a refined node
          endif
!
!  .....if an active node
        endif
!
!  ...end of loop through nodes
      enddo


      call reset_visit

!
! ...second loop through the elements to enforce min rule for the edges wrt faces
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_nodes(mdle, nodesl,norientl)
        call element_order(mdle,norientl, norder)
        etype = NODES(mdle)%type
        nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!  .....collect the order for nodes determined so far in the ELEMENT coordinates
        do j=1,nre+nrf
          nod = nodesl(nrv+j)
          NODES(nod)%visit = norder(nrv+j)
!
!  .......pick up the order determined in the first loop using the element min rule
!
!  .......determine the face order in element coordinates
          select case(NODES(nod)%type)
          case('mdlq')
            norder(j) = NODES(nod)%order
            call decode(norder(j), nordh,nordv)
            select case(norientl(nrv+j))
            case(1,3,4,6); norder(j) = nordv*10+nordh
            end select
          end select
        enddo
!
!    ...modify the order of edges according to the selected rule
!    ...loop through faces
        do jf=1,nrf
!
!   ......determine element edge numbers for the face
          call face_to_edge(etype,jf, ne1,ne2,ne3,ne4)
          nod = nodesl(nrv+nre+jf)
          select case(face_type(etype,jf))
          case('tria')
            norder(ne1) = min(norder(ne1),norder(nre+jf))
            norder(ne2) = min(norder(ne2),norder(nre+jf))
            norder(ne3) = min(norder(ne3),norder(nre+jf))
          case('rect')
            call decode(norder(nre+jf), nordh,nordv)
            norder(ne1) = min(norder(ne1),nordh)
            norder(ne2) = min(norder(ne2),nordv)
            norder(ne3) = min(norder(ne3),nordh)
            norder(ne4) = min(norder(ne4),nordv)
          end select
        enddo
!
!  .....loop through edges
        do je=1,nre
          nod = nodesl(nrv+je)
          if (nod .eq. 464) then 
            write(*,*) '1:nod, NODES(nod)%visit = ', nod, norder(j)
          endif 
          NODES(nod)%visit = min(NODES(nod)%visit,norder(je))
        enddo
!  ...end of loop through the elements
      enddo        


!  ...loop through nodes
      do nod=1,NRNODS

        if (nod .eq. 464) then 
          write(*,*) '2:nod, NODES(nod)%visit = ', nod, NODES(nod)%visit
        endif  

        if(NODES(nod)%type .ne. 'medg') cycle
!
!  .....active node
        if (NODES(nod)%act.eq.1) then
          nord = NODES(nod)%visit
          if (nord.eq.0) cycle
!
!  .......refined node, check order of its sons
          if (NODES(nod)%ref_kind.ne.0) then
!
            do is=1,2
              nods = NODES(nod)%sons(is)
              nord = min(nord,NODES(nods)%visit)
            enddo
!
!...........modify the node order
            if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
! ..........communicate the new order to the (inactive) sons
            do is=1,2
              nods = NODES(nod)%sons(is)
              NODES(nods)%order = nord
            enddo
!
!  .......unrefined active node
          else
!
!    .......modify the node order
            nord = NODES(nod)%visit
            if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
! .......if a refined node
          endif
!
!  .....if an active node
        endif
!
!  ...end of loop through nodes
      enddo
!
!  ...reset visitation flags
      call reset_visit
!
 

!
      end subroutine enforce_min_rule
