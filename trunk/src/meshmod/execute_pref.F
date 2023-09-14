c---------------------------------------------------------------------
c
c   routine name       - execute_pref
c
c---------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine executes p-refinement for a group
c                        of elements
c                        REMARK: this routine executes ONLY isotropic
c                                refinements
c
c   arguments :
c
c     in:
c             List     - a list of elements (active or not)
c             Nlist    - length of the list
c
c-----------------------------------------------------------------------
c
      subroutine execute_pref(List,Nlist)
c
      use data_structure3D
      use element_data
      use constrained_nodes
c
      implicit none
c
      integer, intent(in) :: Nlist
      integer, intent(in) :: List(Nlist)
c
c  ...element nodes
      integer :: nodesl(27), norientl(27)
c
c  ...order for element nodes implied by the order of the middle node
      integer :: norder(19)
c
c  ...miscellanea
      integer :: icase,iel,ie,ip,is,j,mdle,nc,nrs
      integer :: nod,newp,nodp,nods,nrv,nre,nrf,ntype
c
#if DEBUG_MODE
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
        write(*,7010) List(1:Nlist)
 7010   format('execute_pref: List = ',20(/,10i6))
        call pause
      endif
#endif
c
      call reset_visit
c
c  ...Step 1: loop over elements from the list and record
c     the new order of approximation for their nodes
      do iel=1,Nlist
        mdle = List(iel)
c
c   ....skip if the element has been h-refined
        if (NODES(mdle)%ref_kind.ne.0) cycle
c
c  .....determine nodes for the element (active and constrained)
c       and build the data base for the constrained nodes
c       (module constrained_nodes)
        call get_connect_info(mdle, nodesl,norientl)
        ntype = NODES(mdle)%ntype
        nrv = NVERT(ntype); nre = NEDGE(ntype); nrf = NFACE(ntype)
c
c  .....loop through element higher order nodes
C         do j=nrv+1,nrv+nre+nrf+1
!     ...ONLY THROUGH FACES AND THE MIDDLE NODE
        do j=nrv+nre+1,nrv+nre+nrf+1
          nod = Nodesl(j)
c
c  .......if nod is active
          if (Is_active(nod)) then
c
            NODES(nod)%visit = 1
            nrs = NODES(nod)%nr_sons
            do is=1,nrs
              nods = Son(nod,is)
              NODES(nods)%visit = 1
            enddo
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
              NODES(nodp)%visit = 1
c
c  .........mid-face node constrained by an h4-refined face...............
            case(21,22,23,24)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  .........horizontal mid-edge node constrained by an h4-refined face....
            case(26,28)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  ...........parent mid-edge nodes (south,north)
              do ip=1,3,2
                nodp = abs(NFACE_CONS(ip,nc))
                NODES(nodp)%visit = 1
              enddo
c
c  .........vertical mid-edge node constrained by an h4-refined face......
            case(25,27)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  ...........parent mid-edge nodes (east,west)
              do ip=2,4,2
                nodp = abs(NFACE_CONS(ip,nc))
                NODES(nodp)%visit = 1
              enddo
c
c  .........mid-face node constrained by a horizontally h2-refined face...
            case(31,32, 34,35, 61,62)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  .........horizontal mid-edge node constrained by a horizontally
c           h2-refined face...............................................
            case(33,36,63)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  ...........parent mid-edge nodes (south,north)
              do ip=1,3,2
                nodp = abs(NFACE_CONS(ip,nc))
                NODES(nodp)%visit = 1
                if (Is_inactive(nodp)) then
                  nodp = NODES(nodp)%father
                  NODES(nodp)%visit = 1
                endif
              enddo
c
c  .........mid-face node constrained by a vertically h2-refined face.....
            case(41,42, 44,45, 51,52)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  .........vertical mid-edge node constrained by a vertically h2-refined
c           face............................................................
            case(43,46,53)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  ...........parent mid-edge nodes (east,west)
              do ip=2,4,2
                nodp = abs(NFACE_CONS(ip,nc))
                NODES(nodp)%visit = 1
                if (Is_inactive(nodp)) then
                  nodp = NODES(nodp)%father
                  NODES(nodp)%visit = 1
                endif
              enddo
c
c  .........mdlt or medg node constrained by a face.....
            case(71,72,73,74,75,76,77)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  ...........parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
                NODES(nodp)%visit = 1
              enddo
c
c  .........mdlq node constrained by an h2-refined triangular face.....
            case(82,83,84)
c
c  ...........parent mid-face node
              nodp = NFACEC(nc)
              NODES(nodp)%visit = 1
c
c  ...........parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
                NODES(nodp)%visit = 1
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
c  ...Step 2: loop through all elements
      do iel=1,NRELES
        mdle = ELEM_ORDER(iel)
        call elem_nodes(mdle, nodesl, norientl)
c
c  .....determine order of element nodes implied by the middle
c       node order
        call element_order(mdle,norientl, norder)
c
c  .....determine whether the element will be p-refined
        if (NODES(mdle)%visit.eq.0) then
          ntype = NODES(mdle)%ntype
          nrv = NVERT(ntype); nre = NEDGE(ntype); nrf = NFACE(ntype)
c
c  .......loop through element higher order nodes
          do j=1,nre+nrf+1
            nod = Nodesl(nrv+j)
c
c  .........if nod is active
            if (NODES(nod)%visit.eq.1) then
              call find_new_order(nod, newp)
              if (newp.gt.norder(j)) then
                NODES(mdle)%visit=1
                exit
              endif
            endif
          enddo
        endif
c
        if (NODES(mdle)%visit.eq.1) then
          call find_new_order(mdle, newp)
#if DEBUG_MODE
          if (iprint.eq.1) then
            write(*,7050) mdle,newp
 7050       format('execute_pref: REFINING mdle = ',i6,' new p = ',i3)
          endif
#endif
          call nodmod(mdle, newp)
        endif
      enddo

      call reset_visit
c
      end subroutine execute_pref
c
c---------------------------------------------------------------------
c
c   routine name       - find_new_order
c
c---------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine determines the increased order
c                        for a node
c   arguments :
c
c     in:
c             Nod      - node number
c     out:
c             Newp     - new order of approximation
c
c-----------------------------------------------------------------------
c
      subroutine find_new_order(Nod, Newp)
c
      use data_structure3D
      implicit none
c
      integer, intent(in)  :: Nod
      integer, intent(out) :: Newp
c
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
c
      end subroutine find_new_order

