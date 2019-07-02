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
!  .........edge node
            case('medg')
              do is=1,2
                nods = NODES(nod)%sons(is)
                nord = min(nord,NODES(nods)%visit)
!
                if (NODES(nods)%ref_kind.ne.0) then
                  write(*,*) 'enforce_min_rule: INCONSISTENCY 1'
                  stop 1
                endif
              enddo
!
!  ...........modify the node order
              if (nord.ne.NODES(Nod)%order) call nodmod(nod, nord)
!
!  ...........communicate the new order to the (inactive) sons
              do is=1,2
                nods = NODES(nod)%sons(is)
                if (NODES(nods)%act.ne.0) then
                  write(*,*) 'enforce_min_rule: INCONSISTENCY 2'
                  stop 1
                endif
                NODES(nods)%order = nord
              enddo
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
            if (iprint.eq.1) &
            write(*,*) 'enforce_min_rule: CALLING nodmod, nod,nord = ',nod,nord
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
!
!  ...reset visitation flags
      call reset_visit
!
!
      end subroutine enforce_min_rule
