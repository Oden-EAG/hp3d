c----------------------------------------------------------------------
c
c   routine name       - logic_nodes
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine establishes list of nodes
c                        for a modified element
c
c   arguments :
c     in:
c            Mdle      - an element number, same as the middle node
c                        number
c            Nodesl    - (local) element nodes as returned by
c                        get_connect_info
c     out:
c            Nodm      - actual (unconstrained) nodes returned in the
c                        standard order: vertex, mid-edge, mid-face
c                        and middle nodes
c            Nrnodm    - number of nodes for the modified element
c
c----------------------------------------------------------------------
c
      subroutine logic_nodes(Mdle,Nodesl, Nodm,Nrnodm)
c
      use element_data
      use data_structure3D
      use constrained_nodes
c
      implicit none
c
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: Nodesl(27)
      integer, intent(out) :: Nodm(MAXNODM)
      integer, intent(out) :: Nrnodm
c
c  ...element type
      integer :: ntype
c
c  ...local lists of nodes
      integer :: list_v(16),list_e(24),list_f(12)
c
c  ...miscellanea
      integer :: ie,ip,iv,jv,icv,ice,icf,icase,j
      integer :: nvoid,nrnodl,nod,nodp,nc,nce
c
#if DEBUG_MODE
      integer :: iprint = 0
#endif
c
c----------------------------------------------------------------------
c
c  ...initialize output
      Nodm(1:MAXNODM) = 0
c
      ntype = NODES(Mdle)%ntype
c
c  ...number of (local) nodes for the element (- middle node)
      nrnodl = Nvert(ntype)+Nedge(ntype)+Nface(ntype)
c
c  ...establish lists of modified element vertex, edge, and face
c     nodes
      icv=0; ice=0; icf=0
      do j=1,nrnodl
        nod = Nodesl(j)
c
c  .....if nod is active
        if ((Is_active(nod)).or.
     .      (Is_inactive(nod) .and. NODES(nod)%ref_kind.ne.0)) then
c
          if (j.le.Nvert(ntype)) then
            call add_to_list(list_v,16,icv,nod)
          elseif (j.le.Nvert(ntype)+Nedge(ntype)) then
            call add_to_list(list_e,24,ice,nod)
          else
            call add_to_list(list_f,12,icf,nod)
          endif
c
c  .....inactive, i.e. constrained node
        else
c
c  .......identify the constraint case
          call decode2(NODES_CONSTR(j), nc,icase)
c
          select case(icase)
c
c  .......first and second mid-edge node constrained by an edge.....
          case(11,12, 37,38, 47,48)
c
c  .........parent mid-edge node
            nodp = NEDGC(nc)
            call add_to_list(list_e,24,ice,nodp)
c
c  .......vertex node constrained by an edge.....
          case(13,39,49)
c
c  .........parent mid-edge node
            nodp = NEDGC(nc)
            call add_to_list(list_e,24,ice,nodp)
c
c  .........loop through parent vertices
            do iv=1,2
              nodp = NEDG_CONS(iv,nc)
c
c  ...........active (unconstrained) vertex node
              if (nodp.gt.0) then
                call add_to_list(list_v,16,icv,nodp)
c
c  ...........inactive (constrained) vertex node
              else
                call decode2(-nodp, nce,nvoid)
                nodp = NEDGC(nce)
                call add_to_list(list_e,24,ice,nodp)
                do jv=1,2
                  nodp = NEDG_CONS(jv,nce)
                  call add_to_list(list_v,16,icv,nodp)
                enddo
              endif
            enddo
c
c  .......mid-face node constrained by an h4-refined face...............
          case(21,22,23,24)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .......horizontal mid-edge node constrained by an h4-refined face....
          case(26,28)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .........parent mid-edge nodes (south,north)
            do ip=1,3,2
              nodp = iabs(NFACE_CONS(ip,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
c
c  .......vertical mid-edge node constrained by an h4-refined face......
          case(25,27)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .........parent mid-edge nodes (east,west)
            do ip=2,4,2
              nodp = iabs(NFACE_CONS(ip,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
c
c  .......vertex node constrained by an h4-refined face...............
          case(29)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .........parent mid-edge nodes (south,east,north,west)
            do ip=1,4
              nodp = iabs(NFACE_CONS(ip,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
c
c  .........parent vertex dof
            do ip=1,4
              nodp = NFACE_CONS(4+ip,nc)
              call add_to_list(list_v,16,icv,nodp)
            enddo
c
c  .......mid-face node constrained by a horizontally h2-refined face...
          case(31,32, 34,35, 61,62)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .......horizontal mid-edge node constrained by a horizontally
c         h2-refined face...............................................
          case(33,36,63)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .........parent mid-edge nodes (south,north)
            do ip=1,3,2
              nodp = iabs(NFACE_CONS(ip,nc))
C               if (Is_inactive(nodp)) nodp = NODES(nodp)%father
              if (Is_inactive(nodp) .and. NODES(nodp)%ref_kind.eq.0)
     .            nodp = NODES(nodp)%father
              call add_to_list(list_e,24,ice,nodp)
            enddo
c
c  .......mid-face node constrained by a vertically h2-refined face.....
          case(41,42, 44,45, 51,52)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .......vertical mid-edge node constrained by a vertically h2-refined
c         face............................................................
          case(43,46,53)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .........parent mid-edge nodes (east,west)
            do ip=2,4,2
              nodp = iabs(NFACE_CONS(ip,nc))
C               if (Is_inactive(nodp)) nodp = NODES(nodp)%father
              if (Is_inactive(nodp) .and. NODES(nodp)%ref_kind.eq.0)
     .            nodp = NODES(nodp)%father
              call add_to_list(list_e,24,ice,nodp)
            enddo
c
c  .......mdlt or medg node constrained by a face.....
          case(71,72,73,74,75,76,77)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .........parent medg nodes
            do ie=1,3
              nodp = iabs(NFACE_CONS(ie,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
c
c  .......mdlq node constrained by an h2-refined triangular face.....
          case(82,83,84)
c
c  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
c
c  .........parent medg nodes
            do ie=1,3
              nodp = iabs(NFACE_CONS(ie,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
c
          end select
c
c  .....if a constrained node
        endif
c
c  ...end of loop through element nodes
      enddo
c
c  ...save the number of modified element nodes
      Nrnodm = icv+ice+icf+1
      if (Nrnodm.gt.MAXNODM) then
        write(*,7010) Nrnodm,MAXNODM
 7010   format('logic_nodes: Nrnodm,MAXNODM = ',2i4)
        stop 1
      endif
c
c  ...put all nodes on the modified element nodes list
      Nodm(1:icv) = list_v(1:icv)
      Nodm(icv+1:icv+ice) = list_e(1:ice)
      Nodm(icv+ice+1:icv+ice+icf) = list_f(1:icf)
      Nodm(Nrnodm) = Mdle
cccc
cccc  ...clean the visitation flag
ccc      do j=1,Nrnodm-1
ccc        NODES(Nodm(j))%visit = 0
ccc      enddo
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7002) Mdle
 7002   format('logic_nodes: Mdle = ',i6,' MODIFIED ELEMENT NODES = ')
        write(*,7003) Nodm(1:Nrnodm)
 7003   format(20i6)
        call pause
      endif
#endif
c
      end subroutine logic_nodes
c
c----------------------------------------------------------------------
c
c   routine name       - add_to_list
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine adds a node but NO LONGER with
c                        visitation flags
c
c   arguments :
c     in:
c            List      - list of nodes
c            Nlist     - length of the list
c            Ic        - counter
c            Nod       - node number
c     out:
c            List,Ic   - modified, if the node has been added
c
c----------------------------------------------------------------------
c
      subroutine add_to_list(List,Nlist,Ic,Nod)
c
      use data_structure3D
c
      implicit none
c
      integer, intent(in)    :: Nlist
      integer, intent(inout) :: List(Nlist)
      integer, intent(inout) :: Ic
      integer, intent(in)    :: Nod
c
      integer :: number
cccc
ccc      if (NODES(Nod)%visit.eq.0) then
ccc        Ic=Ic+1
ccc        if (Ic.gt.Nlist) then
ccc          write(*,7001) Nlist
ccc 7001     format('add_to_list: LENGTH OF List EXCEEDED, cccNlist = ',i3)
ccc          stop 1
ccc        endif
ccc        List(Ic) = Nod
cccc
cccc  .....raise the visitation flag
ccc        NODES(Nod)%visit=1
ccc      endif
      call locate(Nod,List,Ic, number)
c
      if (number.eq.0) then
c  .....add the node to the list
        Ic=Ic+1
        if (Ic.gt.Nlist) then
          write(*,7001) Nlist
 7001     format('add_to_list: LENGTH OF List EXCEEDED, Nlist = ',i3)
          stop 1
        endif
        List(Ic) = Nod
      endif
c
      end subroutine add_to_list
