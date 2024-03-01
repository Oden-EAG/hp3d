!----------------------------------------------------------------------
!
!   routine name       - logic_nodes
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine establishes list of nodes
!                        for a modified element
!
!   arguments :
!     in:
!            Mdle      - an element number, same as the middle node
!                        number
!            Nodesl    - (local) element nodes as returned by
!                        get_connect_info
!     out:
!            Nodm      - actual (unconstrained) nodes returned in the
!                        standard order: vertex, mid-edge, mid-face
!                        and middle nodes
!            Nrnodm    - number of nodes for the modified element
!
!----------------------------------------------------------------------
!
   subroutine logic_nodes(Mdle,Nodesl, Nodm,Nrnodm)
!
      use element_data
      use data_structure3D
      use constrained_nodes
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: Nodesl(27)
      integer, intent(out) :: Nodm(MAXNODM)
      integer, intent(out) :: Nrnodm
!
!  ...element type
      integer :: ntype
!
!  ...local lists of nodes
      integer :: list_v(16),list_e(24),list_f(12)
!
!  ...miscellanea
      integer :: ie,ip,iv,jv,icv,ice,icf,icase,j
      integer :: nvoid,nrnodl,nod,nodp,nc,nce
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!----------------------------------------------------------------------
!
!  ...initialize output
      Nodm(1:MAXNODM) = 0
!
      ntype = NODES(Mdle)%ntype
!
!  ...number of (local) nodes for the element (- middle node)
      nrnodl = Nvert(ntype)+Nedge(ntype)+Nface(ntype)
!
!  ...establish lists of modified element vertex, edge, and face
!     nodes
      icv=0; ice=0; icf=0
      do j=1,nrnodl
        nod = Nodesl(j)
!
!  .....if nod is active
        if ((Is_active(nod)).or. &
            (Is_inactive(nod) .and. NODES(nod)%ref_kind.ne.0)) then
!
          if (j.le.Nvert(ntype)) then
            call add_to_list(list_v,16,icv,nod)
          elseif (j.le.Nvert(ntype)+Nedge(ntype)) then
            call add_to_list(list_e,24,ice,nod)
          else
            call add_to_list(list_f,12,icf,nod)
          endif
!
!  .....inactive, i.e. constrained node
        else
!
!  .......identify the constraint case
          call decode2(NODES_CONSTR(j), nc,icase)
!
          select case(icase)
!
!  .......first and second mid-edge node constrained by an edge.....
          case(11,12, 37,38, 47,48)
!
!  .........parent mid-edge node
            nodp = NEDGC(nc)
            call add_to_list(list_e,24,ice,nodp)
!
!  .......vertex node constrained by an edge.....
          case(13,39,49)
!
!  .........parent mid-edge node
            nodp = NEDGC(nc)
            call add_to_list(list_e,24,ice,nodp)
!
!  .........loop through parent vertices
            do iv=1,2
              nodp = NEDG_CONS(iv,nc)
!
!  ...........active (unconstrained) vertex node
              if (nodp.gt.0) then
                call add_to_list(list_v,16,icv,nodp)
!
!  ...........inactive (constrained) vertex node
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
!
!  .......mid-face node constrained by an h4-refined face...............
          case(21,22,23,24)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .......horizontal mid-edge node constrained by an h4-refined face....
          case(26,28)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .........parent mid-edge nodes (south,north)
            do ip=1,3,2
              nodp = abs(NFACE_CONS(ip,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
!
!  .......vertical mid-edge node constrained by an h4-refined face......
          case(25,27)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .........parent mid-edge nodes (east,west)
            do ip=2,4,2
              nodp = abs(NFACE_CONS(ip,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
!
!  .......vertex node constrained by an h4-refined face...............
          case(29)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .........parent mid-edge nodes (south,east,north,west)
            do ip=1,4
              nodp = abs(NFACE_CONS(ip,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
!
!  .........parent vertex dof
            do ip=1,4
              nodp = NFACE_CONS(4+ip,nc)
              call add_to_list(list_v,16,icv,nodp)
            enddo
!
!  .......mid-face node constrained by a horizontally h2-refined face...
          case(31,32, 34,35, 61,62)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .......horizontal mid-edge node constrained by a horizontally
!         h2-refined face...............................................
          case(33,36,63)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .........parent mid-edge nodes (south,north)
            do ip=1,3,2
              nodp = abs(NFACE_CONS(ip,nc))
!               if (Is_inactive(nodp)) nodp = NODES(nodp)%father
              if (Is_inactive(nodp) .and. NODES(nodp)%ref_kind.eq.0) &
                  nodp = NODES(nodp)%father
              call add_to_list(list_e,24,ice,nodp)
            enddo
!
!  .......mid-face node constrained by a vertically h2-refined face.....
          case(41,42, 44,45, 51,52)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .......vertical mid-edge node constrained by a vertically h2-refined
!         face............................................................
          case(43,46,53)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .........parent mid-edge nodes (east,west)
            do ip=2,4,2
              nodp = abs(NFACE_CONS(ip,nc))
!               if (Is_inactive(nodp)) nodp = NODES(nodp)%father
              if (Is_inactive(nodp) .and. NODES(nodp)%ref_kind.eq.0) &
                  nodp = NODES(nodp)%father
              call add_to_list(list_e,24,ice,nodp)
            enddo
!
!  .......mdlt or medg node constrained by a face.....
          case(71,72,73,74,75,76,77)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .........parent medg nodes
            do ie=1,3
              nodp = abs(NFACE_CONS(ie,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
!
!  .......mdlq node constrained by an h2-refined triangular face.....
          case(82,83,84)
!
!  .........parent mid-face node
            nodp = NFACEC(nc)
            call add_to_list(list_f,12,icf,nodp)
!
!  .........parent medg nodes
            do ie=1,3
              nodp = abs(NFACE_CONS(ie,nc))
              call add_to_list(list_e,24,ice,nodp)
            enddo
!
          end select
!
!  .....if a constrained node
        endif
!
!  ...end of loop through element nodes
      enddo
!
!  ...save the number of modified element nodes
      Nrnodm = icv+ice+icf+1
      if (Nrnodm.gt.MAXNODM) then
        write(*,7010) Nrnodm,MAXNODM
 7010   format('logic_nodes: Nrnodm,MAXNODM = ',2i4)
        stop 1
      endif
!
!  ...put all nodes on the modified element nodes list
      Nodm(1:icv) = list_v(1:icv)
      Nodm(icv+1:icv+ice) = list_e(1:ice)
      Nodm(icv+ice+1:icv+ice+icf) = list_f(1:icf)
      Nodm(Nrnodm) = Mdle
!!!!
!!!!  ...clean the visitation flag
!!!      do j=1,Nrnodm-1
!!!        NODES(Nodm(j))%visit = 0
!!!      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7002) Mdle
 7002   format('logic_nodes: Mdle = ',i6,' MODIFIED ELEMENT NODES = ')
        write(*,7003) Nodm(1:Nrnodm)
 7003   format(20i6)
        call pause
      endif
#endif
!
   end subroutine logic_nodes
!
!----------------------------------------------------------------------
!
!   routine name       - add_to_list
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine adds a node but NO LONGER with
!                        visitation flags
!
!   arguments :
!     in:
!            List      - list of nodes
!            Nlist     - length of the list
!            Ic        - counter
!            Nod       - node number
!     out:
!            List,Ic   - modified, if the node has been added
!
!----------------------------------------------------------------------
!
   subroutine add_to_list(List,Nlist,Ic,Nod)
!
      use data_structure3D
!
      implicit none
!
      integer, intent(in)    :: Nlist
      integer, intent(inout) :: List(Nlist)
      integer, intent(inout) :: Ic
      integer, intent(in)    :: Nod
!
      integer :: number
!!!!
!!!      if (NODES(Nod)%visit.eq.0) then
!!!        Ic=Ic+1
!!!        if (Ic.gt.Nlist) then
!!!          write(*,7001) Nlist
!!! 7001     format('add_to_list: LENGTH OF List EXCEEDED, &
!!!                               Nlist = ',i3)
!!!          stop 1
!!!        endif
!!!        List(Ic) = Nod
!!!!
!!!!  .....raise the visitation flag
!!!        NODES(Nod)%visit=1
!!!      endif
      call locate(Nod,List,Ic, number)
!
      if (number.eq.0) then
!  .....add the node to the list
        Ic=Ic+1
        if (Ic.gt.Nlist) then
          write(*,7001) Nlist
 7001     format('add_to_list: LENGTH OF List EXCEEDED, Nlist = ',i3)
          stop 1
        endif
        List(Ic) = Nod
      endif
!
   end subroutine add_to_list
