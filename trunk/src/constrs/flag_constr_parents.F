c----------------------------------------------------------------------
c
c   routine name       - flag_constr_parents
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - routine flags ancestors of parents of
c                        constrained nodes of an element that
c                        must be unconstrained in constrained
c                        approximation
c
c   arguments :
c     in:
c            Mdle      - an element number, same as the middle node
c                        number
c
c----------------------------------------------------------------------
c
      subroutine flag_constr_parents(Mdle)
c
      use element_data
      use data_structure3D
      use constrained_nodes
#include "syscom.blk"
c
c  ...element type
      character(len=4) :: type
c
c----------------------------------------------------------------------
c
      iprint = 0
c
      type = NODES(Mdle)%type
c
c  ...number of (local) nodes for the element (- middle node)
      nrnodl = Nvert(type)+Nedge(type)+Nface(type)
c
c  ...loop through element constrained nodes
      do j=1,nrnodl
        if (NODES_CONSTR(j).eq.0) cycle
c
c  .....identify the constraint case
        call decode2(NODES_CONSTR(j), nc,icase)
c
        select case(icase)
c
c  .....first and second mid-edge node constrained by an edge.....
        case(11,12, 37,38, 47,48)
c
c  .......parent mid-edge node
          nodp = NEDGC(nc)
          if (Is_inactive(nodp)) call flag(nodp)
c
c  .....vertex node constrained by an edge.....
        case(13,39,49)
c
c  .......parent mid-edge node
          nodp = NEDGC(nc)
          if (Is_inactive(nodp)) call flag(nodp)
c
c  .......loop through parent vertices
          do iv=1,2
            nodp = NEDG_CONS(iv,nc)
c
c  .........vertex node that is supposed to be active
            if (nodp.gt.0) then
              if (Is_inactive(nodp)) call flag(nodp)
c
c  .........inactive vertex node, just address in NEDGC
            else
              call decode2(-nodp, nce,nvoid)
              nodp = NEDGC(nce)
              if (Is_inactive(nodp)) call flag(nodp)
              do jv=1,2
                nodp = NEDG_CONS(jv,nce)
                if (Is_inactive(nodp)) call flag(nodp)
              enddo
            endif
          enddo
c
c  .....mid-face node constrained by an h4-refined face...............
        case(21,22,23,24)
c
c  .......parent mid-face node
          nodp = NFACEC(nc)
c
c  .......check: parent mid-face node MUST be active
          if (Is_inactive(nodp)) then
            write(*,8001) 1; stop 1
 8001       format('flag_constr_parents: INCONSISTENCY ',i2)
          endif
c
c  .....horizontal mid-edge node constrained by an h4-refined face....
        case(26,28)
c
c  .......check: parent mid-face node MUST be active
          nodp = NFACEC(nc)
          if (Is_inactive(nodp)) then
            write(*,8001) 2; stop 1
          endif
c
c  .......parent mid-edge nodes (south,north)
          do ip=1,3,2
            nodp = iabs(NFACE_CONS(ip,nc))
            if (Is_inactive(nodp)) call flag(nodp)
          enddo
c
c  .....vertical mid-edge node constrained by an h4-refined face......
        case(25,27)
c
c  .......check: parent mid-face node MUST be active
          nodp = NFACEC(nc)
          if (Is_inactive(nodp)) then
            write(*,8001) 2; stop 1
          endif
c
c  .......parent mid-edge nodes (east,west)
          do ip=2,4,2
            nodp = iabs(NFACE_CONS(ip,nc))
            if (Is_inactive(nodp)) call flag(nodp)
          enddo
c
c  .....vertex node constrained by an h4-refined face...............
        case(29)
c
c  .......check: parent mid-face node MUST be active
          nodp = NFACEC(nc)
          if (Is_inactive(nodp)) then
            write(*,8001) 2; stop 1
          endif
c
c  .......parent mid-edge nodes (south,east,north,west)
          do ip=1,4
            nodp = iabs(NFACE_CONS(ip,nc))
            if (Is_inactive(nodp)) call flag(nodp)
          enddo
c
c  .......parent vertex dof
          do ip=1,4
            nodp = NFACE_CONS(4+ip,nc)
            if (Is_inactive(nodp)) call flag(nodp)
          enddo
c
c  .....mid-face node constrained by a horizontally h2-refined face...
        case(31,32, 34,35, 61,62)
c
c  .......parent mid-face node
          nodp = NFACEC(nc)
          if (Is_inactive(nodp)) then
            write(*,8001) 3; stop 1
          endif
c
c  .....horizontal mid-edge node constrained by a horizontally
c       h2-refined face...............................................
        case(33,36,63)
c
c  .......parent mid-face node
          nodp = NFACEC(nc)
          if (Is_inactive(nodp)) then
            write(*,8001) 4; stop 1
          endif
c
c  .......parent mid-edge nodes (south,north)
          do ip=1,3,2
            nodp = iabs(NFACE_CONS(ip,nc))
            if (Is_inactive(nodp)) then
c
c  ...........a constraining edge is allowed
              if (NODES(nodp)%type.eq.'medg') nodp = NODES(nodp)%father
            endif
ccc            if (Is_inactive(nodp)) nodp = NODES(nodp)%father
cccc
cccc  .........check: father MUST be a medg node
ccc            if (NODES(nodp)%type.ne.'medg') then
ccc              write(*,8001) 5; stop 1
ccc            endif
            if (Is_inactive(nodp)) call flag(nodp)
          enddo
c
c  .....mid-face node constrained by a vertically h2-refined face.....
        case(41,42, 44,45, 51,52)
c
c  .......parent mid-face node
          nodp = NFACEC(nc)
          if (Is_inactive(nodp)) then
            write(*,8001) 6; stop 1
          endif
c
c  .....vertical mid-edge node constrained by a vertically h2-refined
c       face............................................................
        case(43,46,53)
c
c  .......parent mid-face node
          nodp = NFACEC(nc)
          if (Is_inactive(nodp)) then
            write(*,8001) 7; stop 1
          endif
c
c  .......parent mid-edge nodes (east,west)
          do ip=2,4,2
            nodp = iabs(NFACE_CONS(ip,nc))
            if (Is_inactive(nodp)) then
c
c  ...........a constraining edge is allowed
              if (NODES(nodp)%type.eq.'medg') nodp = NODES(nodp)%father
            endif
c
ccc  .........check: father MUST be a medg node
ccc            if (NODES(nodp)%type.ne.'medg') then
ccc              write(*,*) 'Mdle,NODES(nodp)%type = ',
ccc     .                    Mdle,NODES(nodp)%type
ccc              write(*,8001) 8
ccc            endif
            if (Is_inactive(nodp)) call flag(nodp)
          enddo
c
c  .....mdlt,medg or mdlq node constrained by a triangular face.....
        case(71,72,73,74,75,76,77,82,83,84)
c
c  .......parent mid-face node
          nodp = NFACEC(nc)
          if (Is_inactive(nodp)) then
            write(*,8001) 9; stop 1
          endif
c
c  .......parent medg nodes
          do ie=1,3
            nodp = iabs(NFACE_CONS(ie,nc))
            if (Is_inactive(nodp)) call flag(nodp)
          enddo
c
        end select
c
c  ...end of loop through constrained nodes
      enddo
c
c
      end subroutine flag_constr_parents
c
c----------------------------------------------------------------------
c
c   routine name       - flag
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - routine flags all ancestors of a node
c
c   arguments :
c     in:
c            Nod       - node number
c
c----------------------------------------------------------------------
c
      subroutine flag(Nod)
c
      use element_data
      use data_structure3D
      use constrained_nodes
#include "syscom.blk"
c
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) Nod
 7001   format('flag: Nod = ',i7)
        call pause
      endif
c
      nfath = NODES(nod)%father
      do while(nfath.gt.0)
        NODES(nfath)%visit=1
        nfath = NODES(nfath)%father
      enddo
c
c
      end subroutine

