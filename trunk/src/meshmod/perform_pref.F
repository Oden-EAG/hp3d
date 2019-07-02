c---------------------------------------------------------------------
c
c   routine name       - perform_pref
c
c---------------------------------------------------------------------
c
c   latest revision    - Aug 17
c
c   purpose            - routine executes p-refinements for a group
c                        of element middle nodes, along with the
c                        corresponding p-refinements of edge and face 
c                        nodes
c
c   arguments :
c
c     in:
c            List_elem - a list of active elements
c            List_ord  - the corresponding orders for the middle
c                        nodes
c            Nlist     - length of both lists
c
c-----------------------------------------------------------------------
c
      subroutine perform_pref(List_elem,List_nord,Nlist)
c
      use data_structure3D
      use element_data
      use constrained_nodes
#include "syscom.blk"
c
      dimension List_elem(Nlist),List_nord(Nlist)
c
      character(len=4) :: etype
c
c  ...element nodes
      dimension nodesl(27), norientl(27)
c
c  ...order for element nodes 
      dimension norder(19)
c
      iprint=0
      if (iprint.eq.1) then
        write(*,7010) List_elem(1:Nlist)
 7010   format('perform_pref: List_elem = ',20(/,10i6))
        write(*,7020) List_nord(1:Nlist)
 7020   format('perform_pref: List_nord = ',20(/,10i6))
        call pause
      endif
c
c  ...Step 1: loop over elements from the list and perform the requested
c     p-refinements
      do iel=1,Nlist
        mdle = List_elem(iel)
        nord = List_nord(iel)
c
c   ....skip if not an active element
        if (NODES(mdle)%ref_kind.ne.0) cycle
c
c  .....enforce the new order for the element middle node
        call nodmod(mdle,nord)

c  .....determine nodes for the element (active and constrained)
c       and build the data base for the constrained nodes
c       (module constrained_nodes)
        call get_connect_info(mdle, nodesl, norientl)
c
c  .....determine the order for element nodes implied by the middle node
        call element_order(mdle,norientl,norder)
c
        etype = NODES(mdle)%type
        nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
c
c
c  .....loop through element higher order nodes
        do j=nrv+1,nrv+nre+nrf
c
          nod = Nodesl(j)
c
c  .......new order for the node 
          newp = norder(j-nrv) 
C          
          if (NODES(nod)%type.eq.'mdlq') call decode(newp, newph,newpv)
c
c  .......if nod is active
          if (NODES(nod)%act.eq.1) then
c
c  .........raise the flag to enforce the order
            NODES(nod)%visit = 1
            call enforce_order(nod, newp)
c
c  .......inactive, i.e. constrained node
          else
c
c  .........identify the constraint case
            call decode2(NODES_CONSTR(j), nc,icase)
c
            select case(icase)
c
c  .........first and second mid-edge node constrained by an edge.....
            case(11,12, 37,38, 47,48)
c
c  ...........parent mid-edge node
              nodp = NEDGC(nc)
              call enforce_order(nodp, newp)
c
c  .........mid-face node constrained by an h4-refined face...............
            case(21,22,23,24)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp)
c
c  .........horizontal mid-edge node constrained by an h4-refined face....
            case(26,28)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call decode(NODES(nodp)%order, nordh,nordv)
              call enforce_order(nodp, newp*10+nordv)
c
c  ...........parent mid-edge nodes (south,north)
              do ip=1,3,2
                nodp = iabs(NFACE_CONS(ip,nc))
                call enforce_order(nodp, newp)
              enddo
c
c  .........vertical mid-edge node constrained by an h4-refined face......
            case(25,27)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call decode(NODES(nodp)%order, nordh,nordv)
              call enforce_order(nodp, nordh*10+newp)
c
c  ...........parent mid-edge nodes (east,west)
              do ip=2,4,2
                nodp = iabs(NFACE_CONS(ip,nc))
                call enforce_order(nodp, newp)
              enddo
c
c  .........mid-face node constrained by a horizontally h2-refined face...
            case(31,32, 34,35, 61,62)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp)
c
c  .........horizontal mid-edge node constrained by a horizontally
c           h2-refined face...............................................
            case(33,36,63)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call decode(NODES(nodp)%order, nordh,nordv)
              call enforce_order(nodp, newp*10+nordv)
c
c  ...........parent mid-edge nodes (south,north)
              do ip=1,3,2
                nodp = iabs(NFACE_CONS(ip,nc))
                call enforce_order(nodp, newp)
                if (NODES(nodp)%act.eq.0) then
                  nodp = NODES(nodp)%father
                  call enforce_order(nodp, newp)
                endif
              enddo
c
c  .........mid-face node constrained by a vertically h2-refined face.....
            case(41,42, 44,45, 51,52)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp)
c
c  .........vertical mid-edge node constrained by a vertically h2-refined
c           face............................................................
            case(43,46,53)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call decode(NODES(nodp)%order, nordh,nordv)
              call enforce_order(nodp, nordh*10+newp)
c
c  ...........parent mid-edge nodes (east,west)
              do ip=2,4,2
                nodp = iabs(NFACE_CONS(ip,nc))
                call enforce_order(nodp, newp)
                if (NODES(nodp)%act.eq.0) then
                  nodp = NODES(nodp)%father
                  call enforce_order(nodp, newp)
                endif
              enddo
c
c  .........mdlt or medg node constrained by a face.....
            case(71,72,73,74,75,76,77)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp)
c
c  ...........parent medg nodes
              do ie=1,3
                nodp = iabs(NFACE_CONS(ie,nc))
                call enforce_order(nodp, newp)
              enddo
c
c  .........mdlq node constrained by an h2-refined triangular face.....
            case(82,83,84)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp*10+newp)
c
c  ...........parent medg nodes
              do ie=1,3
                nodp = iabs(NFACE_CONS(ie,nc))
               call enforce_order(nodp, newp)
              enddo
c
            end select
c
c  .......if a constrained node
          endif
c
c  .....end of loop through element nodes
        enddo
c
c  ...end of loop through elements
      enddo
c
c-----------------------------------------------------------------------
c
c  ...Step 2: loop through all active elements
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_nodes(mdle, nodesl, norientl)
c
c  .....determine order of element nodes implied by the middle
c       node order
        call element_order(mdle,norientl, norder)
c
c  .....determine whether the element will be p-refined
        if (NODES(mdle)%visit.eq.0) then
          etype = NODES(mdle)%type
          nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
c
c  .......loop through element higher order nodes
          do j=1,nre+nrf+1
            nod = Nodesl(nrv+j)
c
c  .........if the node order is protected
            if (NODES(nod)%visit.eq.1) then
c
c  ...........trade the order implied by the middle node for the
c             actual node order
              norder(j) = NODES(nod)%order
c
c  ...........indicate that the middle node has to be modified
              NODES(mdle)%visit=1
            endif
          enddo
        endif
c
c  .....the middle node has to be modified
        if (NODES(mdle)%visit.eq.1) then
          call element_middle_node_order(Mdle,norientl,norder, nordm)
          if (iprint.eq.1) then
            write(*,7050) mdle,nordm
 7050       format('perform_pref: REFINING mdle = ',i6,' new p = ',i3)
          endif
          call nodmod(mdle, nordm)
        endif
      enddo
c
c  ...reset the visitation flags
      call reset_visit
c
      end subroutine perform_pref
c
c---------------------------------------------------------------------
c
c   routine name       - enforce_order
c
c---------------------------------------------------------------------
c
c   latest revision    - Aug 17
c
c   purpose            - routine modifies order for a node enforcing
c                        the min rule if the node is flagged as
c                        'protected'
c   arguments :
c
c     in:
c             Nod      - node number
c     out:
c             Newp     - new order of approximation
c
c-----------------------------------------------------------------------
c
      subroutine enforce_order(Nod,Newp)
c
      use data_structure3D
c
      select case(NODES(Nod)%visit)
      case(0)
        nord = Newp
        NODES(Nod)%visit = 1
      case default
        select case(NODES(Nod)%type)
        case('medg','mdlt')
          nord = min(Newp,NODES(Nod)%order)
        case('mdlq')
          call decode(Newp, newph,newpv)
          call decode(NODES(Nod)%order, nordh,nordv)
          nord = min(newph,nordh)*10 + min(newpv,nordv)
        end select
      end select
      call nodmod(Nod,nord)
c
c
      end subroutine enforce_order
