!---------------------------------------------------------------------
!
!   routine name       - perform_pref
!
!---------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine executes p-refinements for a group
!                        of element middle nodes, along with the
!                        corresponding p-refinements of edge and face
!                        nodes
!
!   arguments :
!
!     in:
!            List_elem - a list of active elements
!            List_ord  - the corresponding orders for the middle
!                        nodes
!            Nlist     - length of both lists
!
!    REMARK: routine appears to call nodmod more often than necessary
!            need to check correctness of the routine
!
!-----------------------------------------------------------------------
!
      subroutine perform_pref(List_elem,List_nord,Nlist)
!
      use data_structure3D
      use element_data
      use constrained_nodes
!
      implicit none
!
      integer, intent(in) :: Nlist
      integer, intent(in) :: List_elem(Nlist), List_nord(Nlist)
!
!  ...element nodes
      integer :: nodesl(27), norientl(27)
!
!  ...order for element nodes
      integer :: norder(19)
!
      integer :: ie,iel,ip,j,icase,nc,mdle,ntype,nrv,nre,nrf
      integer :: nod,nodp,newp,newph,newpv,nord,nordh,nordv,nordm
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
!
      if (iprint.eq.1) then
        write(*,7010) List_elem(1:Nlist)
 7010   format('perform_pref: List_elem = ',20(/,10i6))
        write(*,7020) List_nord(1:Nlist)
 7020   format('perform_pref: List_nord = ',20(/,10i6))
        call pause
      endif
#endif
!
!  ...Step 1: loop over elements from the list and perform the requested
!     p-refinements
      do iel=1,Nlist
        mdle = List_elem(iel)
        nord = List_nord(iel)
!
!   ....skip if not an active element
        if (NODES(mdle)%ref_kind.ne.0) cycle
!
!  .....enforce the new order for the element middle node
        call nodmod(mdle,nord)

!  .....determine nodes for the element (active and constrained)
!       and build the data base for the constrained nodes
!       (module constrained_nodes)
        call get_connect_info(mdle, nodesl,norientl)
!
!  .....determine the order for element nodes implied by the middle node
        call element_order(mdle,norientl, norder)
!
        ntype = NODES(mdle)%ntype
        nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!
!  .....loop through element higher order nodes
        do j=nrv+1,nrv+nre+nrf
!
          nod = nodesl(j)
!
!  .......new order for the node
          newp = norder(j-nrv)
!
          if (NODES(nod)%ntype.eq.MDLQ) call decode(newp, newph,newpv)
!
!  .......if nod is active
          if (Is_active(nod)) then
!
!  .........raise the flag to enforce the order
            NODES(nod)%visit = 1
            call enforce_order(nod, newp)
!
!  .......inactive, i.e. constrained node
          else
!
!  .........identify the constraint case
            call decode2(NODES_CONSTR(j), nc,icase)
!
            select case(icase)
!
!  .........first and second mid-edge node constrained by an edge.....
            case(11,12, 37,38, 47,48)
!
!  ...........parent mid-edge node
              nodp = NEDGC(nc)
              call enforce_order(nodp, newp)
!
!  .........mid-face node constrained by an h4-refined face...............
            case(21,22,23,24)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp)
!
!  .........horizontal mid-edge node constrained by an h4-refined face....
            case(26,28)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call decode(NODES(nodp)%order, nordh,nordv)
              call enforce_order(nodp, newp*10+nordv)
!
!  ...........parent mid-edge nodes (south,north)
              do ip=1,3,2
                nodp = abs(NFACE_CONS(ip,nc))
                call enforce_order(nodp, newp)
              enddo
!
!  .........vertical mid-edge node constrained by an h4-refined face......
            case(25,27)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call decode(NODES(nodp)%order, nordh,nordv)
              call enforce_order(nodp, nordh*10+newp)
!
!  ...........parent mid-edge nodes (east,west)
              do ip=2,4,2
                nodp = abs(NFACE_CONS(ip,nc))
                call enforce_order(nodp, newp)
              enddo
!
!  .........mid-face node constrained by a horizontally h2-refined face...
            case(31,32, 34,35, 61,62)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp)
!
!  .........horizontal mid-edge node constrained by a horizontally
!           h2-refined face...............................................
            case(33,36,63)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call decode(NODES(nodp)%order, nordh,nordv)
              call enforce_order(nodp, newp*10+nordv)
!
!  ...........parent mid-edge nodes (south,north)
              do ip=1,3,2
                nodp = abs(NFACE_CONS(ip,nc))
                call enforce_order(nodp, newp)
                if (Is_inactive(nodp)) then
                  nodp = NODES(nodp)%father
                  call enforce_order(nodp, newp)
                endif
              enddo
!
!  .........mid-face node constrained by a vertically h2-refined face.....
            case(41,42, 44,45, 51,52)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp)
!
!  .........vertical mid-edge node constrained by a vertically h2-refined
!           face............................................................
            case(43,46,53)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call decode(NODES(nodp)%order, nordh,nordv)
              call enforce_order(nodp, nordh*10+newp)
!
!  ...........parent mid-edge nodes (east,west)
              do ip=2,4,2
                nodp = abs(NFACE_CONS(ip,nc))
                call enforce_order(nodp, newp)
                if (Is_inactive(nodp)) then
                  nodp = NODES(nodp)%father
                  call enforce_order(nodp, newp)
                endif
              enddo
!
!  .........mdlt or medg node constrained by a face.....
            case(71,72,73,74,75,76,77)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp)
!
!  ...........parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
                call enforce_order(nodp, newp)
              enddo
!
!  .........mdlq node constrained by an h2-refined triangular face.....
            case(82,83,84)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              call enforce_order(nodp, newp*10+newp)
!
!  ...........parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
               call enforce_order(nodp, newp)
              enddo
!
            end select
!
!  .......if a constrained node
          endif
!
!  .....end of loop through element nodes
        enddo
!
!  ...end of loop through elements
      enddo
!
!-----------------------------------------------------------------------
!
!  ...Step 2: loop through all active elements
      do iel=1,NRELES
        mdle = ELEM_ORDER(iel)
        call elem_nodes(mdle, nodesl,norientl)
!
!  .....determine order of element nodes implied by the middle
!       node order
        call element_order(mdle,norientl, norder)
!
!  .....determine whether the element will be p-refined
        if (NODES(mdle)%visit.eq.0) then
          ntype = NODES(mdle)%ntype
          nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!  .......loop through element higher order nodes
          do j=1,nre+nrf+1
            nod = nodesl(nrv+j)
!
!  .........if the node order is protected
            if (NODES(nod)%visit.eq.1) then
!
!  ...........trade the order implied by the middle node for the
!             actual node order
              norder(j) = NODES(nod)%order
!
!  ...........indicate that the middle node has to be modified
              NODES(mdle)%visit=1
            endif
          enddo
        endif
!
!  .....the middle node has to be modified
        if (NODES(mdle)%visit.eq.1) then
          call element_middle_node_order(Mdle,norientl,norder, nordm)
#if HP3D_DEBUG
          if (iprint.eq.1) then
            write(*,7050) mdle,nordm
 7050       format('perform_pref: REFINING mdle = ',i6,' new p = ',i3)
          endif
#endif
          call nodmod(mdle, nordm)
        endif
      enddo
!
!  ...reset the visitation flags
      call reset_visit
!
      end subroutine perform_pref
!
!---------------------------------------------------------------------
!
!   routine name       - enforce_order
!
!---------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine modifies order for a node enforcing
!                        the min rule if the node is flagged as
!                        'protected'
!   arguments :
!
!     in:
!             Nod      - node number
!             Newp     - new order of approximation
!
!-----------------------------------------------------------------------
!
      subroutine enforce_order(Nod,Newp)
!
      use data_structure3D
!
      implicit none
!
      integer, intent(in) :: Nod
      integer, intent(in) :: Newp
!
      integer :: nord,nordh,nordv,newph,newpv
!
      select case(NODES(Nod)%visit)
      case(0)
        nord = Newp
        NODES(Nod)%visit = 1
      case default
        select case(NODES(Nod)%ntype)
        case(MEDG,MDLT)
          nord = min(Newp,NODES(Nod)%order)
        case(MDLQ)
          call decode(Newp, newph,newpv)
          call decode(NODES(Nod)%order, nordh,nordv)
          nord = min(newph,nordh)*10 + min(newpv,nordv)
        end select
      end select
      call nodmod(Nod,nord)
!
      end subroutine enforce_order
