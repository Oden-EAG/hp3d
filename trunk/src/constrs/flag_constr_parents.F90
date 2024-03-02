!----------------------------------------------------------------------
!
!   routine name       - flag_constr_parents
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine flags ancestors of parents of
!                        constrained nodes of an element that
!                        must be unconstrained in constrained
!                        approximation
!
!   arguments :
!     in:
!            Mdle      - an element number, same as the middle node
!                        number
!            Nodesl    - list of element nodes
!
!----------------------------------------------------------------------
!
   subroutine flag_constr_parents(Mdle,Nodesl)
!
      use element_data
      use data_structure3D
      use constrained_nodes
!
      implicit none
!
      integer, intent(in) :: Mdle
      integer, intent(in) :: Nodesl(27)
!
!  ...element type
      integer :: ntype
!
!  ...misc
      integer :: icase,ie,ip,iv,j,jv
      integer :: nrnodl,nod,nodp,nc,nce,nvoid
!
!----------------------------------------------------------------------
!
      ntype = NODES(Mdle)%ntype
!
!  ...number of (local) nodes for the element (- middle node)
      nrnodl = Nvert(ntype)+Nedge(ntype)+Nface(ntype)
!
!  ...loop through element constrained nodes
      do j=1,nrnodl
!
         if (NODES_CONSTR(j).eq.0) then
!        ...check if an active node
            nod = Nodesl(j)
            if (.not.NODES(nod)%act) call flag(nod)
            cycle
         endif
!
!     ...identify the constraint case
         call decode2(NODES_CONSTR(j), nc,icase)
!
         select case(icase)
!
!     ...first and second mid-edge node constrained by an edge
         case(11,12, 37,38, 47,48)
!
!        ...parent mid-edge node
            nodp = NEDGC(nc)
            if (Is_inactive(nodp)) call flag(nodp)
!
!     ...vertex node constrained by an edge
         case(13,39,49)
!
!        ...parent mid-edge node
            nodp = NEDGC(nc)
            if (Is_inactive(nodp)) call flag(nodp)
!
!        ...loop through parent vertices
            do iv=1,2
               nodp = NEDG_CONS(iv,nc)
!
!           ...vertex node that is supposed to be active
               if (nodp.gt.0) then
                  if (Is_inactive(nodp)) call flag(nodp)
!
!           ...inactive vertex node, just address in NEDGC
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
!
!     ...mid-face node constrained by an h4-refined face
         case(21,22,23,24)
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
!
!        ...check: parent mid-face node MUST be active
            if (Is_inactive(nodp)) then
               write(*,8001) 1; stop 1
 8001          format('flag_constr_parents: INCONSISTENCY ',i2)
            endif
!
!     ...horizontal mid-edge node constrained by an h4-refined face
         case(26,28)
!
!        ...check: parent mid-face node MUST be active
            nodp = NFACEC(nc)
            if (Is_inactive(nodp)) then
               write(*,8001) 2; stop 1
            endif
!
!        ...parent mid-edge nodes (south,north)
            do ip=1,3,2
               nodp = iabs(NFACE_CONS(ip,nc))
               if (Is_inactive(nodp)) call flag(nodp)
            enddo
!
!     ...vertical mid-edge node constrained by an h4-refined face
         case(25,27)
!
!        ...check: parent mid-face node MUST be active
            nodp = NFACEC(nc)
            if (Is_inactive(nodp)) then
               write(*,8001) 2; stop 1
            endif
!
!        ...parent mid-edge nodes (east,west)
            do ip=2,4,2
               nodp = iabs(NFACE_CONS(ip,nc))
               if (Is_inactive(nodp)) call flag(nodp)
            enddo
!
!     ...vertex node constrained by an h4-refined face
         case(29)
!
!        ...check: parent mid-face node MUST be active
            nodp = NFACEC(nc)
            if (Is_inactive(nodp)) then
               write(*,8001) 2; stop 1
            endif
!
!        ...parent mid-edge nodes (south,east,north,west)
            do ip=1,4
               nodp = iabs(NFACE_CONS(ip,nc))
               if (Is_inactive(nodp)) call flag(nodp)
            enddo
!
!        ...parent vertex dof
            do ip=1,4
               nodp = NFACE_CONS(4+ip,nc)
               if (Is_inactive(nodp)) call flag(nodp)
            enddo
!
!     ...mid-face node constrained by a horizontally h2-refined face
         case(31,32, 34,35, 61,62)
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            if (Is_inactive(nodp)) then
               write(*,8001) 3; stop 1
            endif
!
!     ...horizontal mid-edge node constrained by a horizontally
!        h2-refined face
         case(33,36,63)
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            if (Is_inactive(nodp)) then
               write(*,8001) 4; stop 1
            endif
!
!        ...parent mid-edge nodes (south,north)
            do ip=1,3,2
               nodp = iabs(NFACE_CONS(ip,nc))
!
               if (Is_inactive(nodp)) then
!              ...a constraining edge is allowed
                  if (NODES(nodp)%ntype.eq.MEDG) nodp = NODES(nodp)%father
               endif
!
               if (Is_inactive(nodp)) call flag(nodp)
            enddo
!
!     ...mid-face node constrained by a vertically h2-refined face
         case(41,42, 44,45, 51,52)
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            if (Is_inactive(nodp)) then
               write(*,8001) 6; stop 1
            endif
!
!     ...vertical mid-edge node constrained by a vertically h2-refined face
         case(43,46,53)
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            if (Is_inactive(nodp)) then
               write(*,8001) 7; stop 1
            endif
!
!        ...parent mid-edge nodes (east,west)
            do ip=2,4,2
               nodp = iabs(NFACE_CONS(ip,nc))
               if (Is_inactive(nodp)) then
!              ...a constraining edge is allowed
                  if (NODES(nodp)%ntype.eq.MEDG) nodp = NODES(nodp)%father
               endif
!
               if (Is_inactive(nodp)) call flag(nodp)
            enddo
!
!     ...mdlt,medg or mdlq node constrained by a triangular face.....
         case(71,72,73,74,75,76,77,82,83,84)
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            if (Is_inactive(nodp)) then
               write(*,8001) 9; stop 1
            endif
!
!        ...parent medg nodes
            do ie=1,3
               nodp = iabs(NFACE_CONS(ie,nc))
               if (Is_inactive(nodp)) call flag(nodp)
            enddo
!
         end select
!
!  ...end of loop through constrained nodes
      enddo
!
!
   end subroutine flag_constr_parents
!
!----------------------------------------------------------------------
!
!   routine name       - flag
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine flags all ancestors of a node
!
!   arguments :
!     in:
!            Nod       - node number
!
!----------------------------------------------------------------------
!
   subroutine flag(Nod)
!
      use data_structure3D,   only : NODES
      use bitvisit
!
      implicit none
!
      integer, intent(in) :: Nod
      integer :: nfath
!
#if HP3D_DEBUG
      integer :: iprint
!
      iprint=0
      if (iprint.eq.1) then
         write(*,7001) Nod
 7001    format('flag: Nod = ',i7)
         call pause
      endif
#endif
!
      nfath = NODES(Nod)%father
      do while(nfath.gt.0)
         call visit(nfath)
         nfath = NODES(nfath)%father
      enddo
!
   end subroutine flag
