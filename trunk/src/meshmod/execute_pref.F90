!---------------------------------------------------------------------
!
!   routine name       - execute_pref
!
!---------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine executes p-refinement for a group
!                        of elements
!                        REMARK: this routine executes ONLY isotropic
!                                refinements
!
!   arguments :
!
!     in:
!             List     - a list of elements (active or not)
!             Nlist    - length of the list
!
!-----------------------------------------------------------------------
!
      subroutine execute_pref(List,Nlist)
!
      use data_structure3D
      use element_data
      use constrained_nodes
!
      implicit none
!
      integer, intent(in) :: Nlist
      integer, intent(in) :: List(Nlist)
!
!  ...element nodes
      integer :: nodesl(27), norientl(27)
!
!  ...order for element nodes implied by the order of the middle node
      integer :: norder(19)
!
!  ...miscellanea
      integer :: icase,iel,ie,ip,is,j,mdle,nc,nrs
      integer :: nod,newp,nodp,nods,nrv,nre,nrf,ntype
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
        write(*,7010) List(1:Nlist)
 7010   format('execute_pref: List = ',20(/,10i6))
        call pause
      endif
#endif
!
      call reset_visit
!
!  ...Step 1: loop over elements from the list and record
!     the new order of approximation for their nodes
      do iel=1,Nlist
        mdle = List(iel)
!
!   ....skip if the element has been h-refined
        if (NODES(mdle)%ref_kind.ne.0) cycle
!
!  .....determine nodes for the element (active and constrained)
!       and build the data base for the constrained nodes
!       (module constrained_nodes)
        call get_connect_info(mdle, nodesl,norientl)
        ntype = NODES(mdle)%ntype
        nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!  .....loop through element higher order nodes
!         do j=nrv+1,nrv+nre+nrf+1
!     ...ONLY THROUGH FACES AND THE MIDDLE NODE
        do j=nrv+nre+1,nrv+nre+nrf+1
          nod = Nodesl(j)
!
!  .......if nod is active
          if (Is_active(nod)) then
!
            NODES(nod)%visit = 1
            nrs = NODES(nod)%nr_sons
            do is=1,nrs
              nods = Son(nod,is)
              NODES(nods)%visit = 1
            enddo
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
              NODES(nodp)%visit = 1
!
!  .........mid-face node constrained by an h4-refined face...............
            case(21,22,23,24)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  .........horizontal mid-edge node constrained by an h4-refined face....
            case(26,28)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  ...........parent mid-edge nodes (south,north)
              do ip=1,3,2
                nodp = abs(NFACE_CONS(ip,nc))
                NODES(nodp)%visit = 1
              enddo
!
!  .........vertical mid-edge node constrained by an h4-refined face......
            case(25,27)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  ...........parent mid-edge nodes (east,west)
              do ip=2,4,2
                nodp = abs(NFACE_CONS(ip,nc))
                NODES(nodp)%visit = 1
              enddo
!
!  .........mid-face node constrained by a horizontally h2-refined face...
            case(31,32, 34,35, 61,62)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  .........horizontal mid-edge node constrained by a horizontally
!           h2-refined face...............................................
            case(33,36,63)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  ...........parent mid-edge nodes (south,north)
              do ip=1,3,2
                nodp = abs(NFACE_CONS(ip,nc))
                NODES(nodp)%visit = 1
                if (Is_inactive(nodp)) then
                  nodp = NODES(nodp)%father
                  NODES(nodp)%visit = 1
                endif
              enddo
!
!  .........mid-face node constrained by a vertically h2-refined face.....
            case(41,42, 44,45, 51,52)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  .........vertical mid-edge node constrained by a vertically h2-refined
!           face............................................................
            case(43,46,53)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  ...........parent mid-edge nodes (east,west)
              do ip=2,4,2
                nodp = abs(NFACE_CONS(ip,nc))
                NODES(nodp)%visit = 1
                if (Is_inactive(nodp)) then
                  nodp = NODES(nodp)%father
                  NODES(nodp)%visit = 1
                endif
              enddo
!
!  .........mdlt or medg node constrained by a face.....
            case(71,72,73,74,75,76,77)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  ...........parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
                NODES(nodp)%visit = 1
              enddo
!
!  .........mdlq node constrained by an h2-refined triangular face.....
            case(82,83,84)
!
!  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
!
!  ...........parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
                NODES(nodp)%visit = 1
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
!  ...Step 2: loop through all elements
      do iel=1,NRELES
        mdle = ELEM_ORDER(iel)
        call elem_nodes(mdle, nodesl, norientl)
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
            nod = Nodesl(nrv+j)
!
!  .........if nod is active
            if (NODES(nod)%visit.eq.1) then
              call find_new_order(nod, newp)
              if (newp.gt.norder(j)) then
                NODES(mdle)%visit=1
                exit
              endif
            endif
          enddo
        endif
!
        if (NODES(mdle)%visit.eq.1) then
          call find_new_order(mdle, newp)
#if HP3D_DEBUG
          if (iprint.eq.1) then
            write(*,7050) mdle,newp
 7050       format('execute_pref: REFINING mdle = ',i6,' new p = ',i3)
          endif
#endif
          call nodmod(mdle, newp)
        endif
      enddo

      call reset_visit
!
      end subroutine execute_pref
!
!---------------------------------------------------------------------
!
!   routine name       - find_new_order
!
!---------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine determines the increased order
!                        for a node
!   arguments :
!
!     in:
!             Nod      - node number
!     out:
!             Newp     - new order of approximation
!
!-----------------------------------------------------------------------
!
      subroutine find_new_order(Nod, Newp)
!
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Nod
      integer, intent(out) :: Newp
!
      select case(NODES(Nod)%ntype)
      case(VERT)
        Newp = 1
      case(MEDG,MDLT,MDLN,MDLD)
        Newp = NODES(Nod)%order+1
      case(MDLQ,MDLP)
        Newp = NODES(Nod)%order+11
      case(MDLB)
        Newp = NODES(Nod)%order+111
      end select
!
      end subroutine find_new_order

