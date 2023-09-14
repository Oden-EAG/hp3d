!----------------------------------------------------------------------
!
!   routine name       - enforce_max_rule
!
!----------------------------------------------------------------------
!> @brief routine enforces the max rule for a FE mesh
!> @remark THIS ROUTINE HAS BEEN OBSERVED TO REFINE (INCREASE P)
!!         UNNECESSARILY WHEN THE MESH HAS ANISOTROPIC ORDER OF APPROX.
!> @date Feb 2023
!----------------------------------------------------------------------
subroutine enforce_max_rule
!
   use data_structure3D
   use refinements
   use constrained_nodes
!
   implicit none
!
!..work space for elem_nodes
   integer :: nodesl(27),norientl(27)
!
!..order for element nodes implied by the order of the middle node
   integer :: norder(19)
!
   integer :: iel,mdle,nre,nrf,nrv,i,j,nod,nord,ntype,nrsons
   integer :: icase,is,jf,nc,ne1,ne2,ne3,ne4,nodp,nods
   integer :: nordh,nordh1,nordhs,nordhl,nordhv
   integer :: nordv,nordv1,nordvs,nordvl
!
!----------------------------------------------------------------------
!
!..reset visitation flags
   call reset_visit
!
!----------------------------------------------------------------------
!        STEP 1: maximum rule for edges and faces wrt elements
!----------------------------------------------------------------------
!
!..loop through elements in the current mesh
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      ntype = NODES(Mdle)%ntype
      nrv = NVERT(ntype); nre = NEDGE(ntype); nrf = NFACE(ntype)
      call get_connect_info(mdle, nodesl,norientl)
      call element_order(mdle,norientl, norder)
!
!  ...loop through edges and faces
      do j=1,nre+nrf
         i=nrv+j
         nod = nodesl(i); nord = norder(j)
!
         call save_max_order(nod,nord)
!
!     ...if a constrained node
         if (Is_inactive(nod)) then
!
!        ...identify the constraint case
            call decode2(NODES_CONSTR(i), nc,icase)
            select case(icase)
!
!        ...edge constrained by an edge
            case(11,12,37,38,47,48)
               nodp = NEDGC(nc)

               call save_max_order(nodp,nord)
!
!        ...horizontal edge constrained by a face
            case(26,28,33,36,63)
               nodp = NFACEC(nc)

               call save_max_order(nodp,nord*10+1)
!
!        ...vertical edge constrained by a face
            case(25,27,43,46,53)
               nodp = NFACEC(nc)
               call save_max_order(nodp,1*10+nord)
!
!        ...face node constrained by a face
            case(21,22,23,24,31,32,34,35,41,42,44,45,51,52,61,62)
               nodp = NFACEC(nc)
               call save_max_order(nodp,nord)
            end select
         endif
      enddo
   enddo
!
!
!----------------------------------------------------------------------
!                 STEP 2: Modify faces (min rule for edges wrt faces)
!----------------------------------------------------------------------
!
!..loop through the elements
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call elem_nodes(mdle, nodesl,norientl)
      ntype = NODES(Mdle)%ntype
      nrv = NVERT(ntype); nre = NEDGE(ntype); nrf = NFACE(ntype)
!
!  ...collect the order for nodes determined so far in the ELEMENT coordinates
      do j=1,nre+nrf
         nod = nodesl(nrv+j)
!
!     ...pick up the order determined in the first loop using the element max rule
         norder(j) = NODES(nod)%visit
!
!     ...determaxe the face order in element coordinates
         select case(NODES(nod)%ntype)
         case(MDLQ)
            call decode(norder(j), nordh,nordv)
            select case(norientl(nrv+j))
            case(1,3,4,6); norder(j) = nordv*10+nordh
            end select
         end select
      enddo
!
!  ...modify the order of the faces according to the edges
!  ...loop through faces
      do jf=1,nrf
!
!     ...determine local edge numbers for the face
         call face_to_edge(ntype,jf, ne1,ne2,ne3,ne4)
         nod = nodesl(nrv+nre+jf)
         select case(face_type(ntype,jf))
         case(TRIA)
            norder(nre+jf) = max(norder(nre+jf),norder(ne1),norder(ne2),norder(ne3))
         case(RECT)
            call decode(norder(nre+jf), nordh,nordv)
            nordh = max(nordh,norder(ne1),norder(ne3))
            nordv = max(nordv,norder(ne2),norder(ne4))
            norder(nre+jf) = nordh*10+nordv
         end select
      enddo
!
!  ...loop through faces
      do jf=1,nrf
         nod = nodesl(nrv+nre+jf)
         select case(face_type(ntype,jf))
         case(TRIA)
            NODES(nod)%visit = max(NODES(nod)%visit,norder(nre+jf))
         case(RECT)
            call decode(NODES(nod)%visit, nordh,nordv)
            call decode(norder(nre+jf), nordhl,nordvl)
            select case(norientl(nrv+nre+jf))
            case(0,2,5,7)
               nordh = max(nordh,nordhl); nordv = max(nordv,nordvl)
            case(1,3,4,6)
               nordh = max(nordh,nordvl); nordv = max(nordv,nordhl)
            end select
            nordhv = max(nordh,nordv)
            NODES(nod)%visit = nordhv*10+nordhv
            ! NODES(nod)%visit = nordh*10+nordv
         end select
      enddo
!
!..end of loop through the elements
   enddo
!
!..loop through nodes
   do nod=1,NRNODS
!
!  ...active node
      if (Is_active(nod)) then
         nord = NODES(nod)%visit
         if (nord.eq.0) cycle
!
!     ...refined node, check order of its sons
         if (NODES(nod)%ref_kind.ne.0) then
            select case(NODES(nod)%ntype)
!
!        ...edge node
            case(MEDG)
               do is=1,2
                  !nods = NODES(nod)%sons(is)
                  nods = Son(nod,is)
                  nord = max(nord,NODES(nods)%visit)
               enddo
!
!        ...modify the node order
            ! write(*,*) 'call nodmod: nod, nord = ', nod, nord
            call nodmod(nod, nord)
!
!        ...communicate the new order to the (inactive) sons
            do is=1,2
               !nods = NODES(nod)%sons(is)
               nods = Son(nod,is)
               NODES(nods)%order = nord
            enddo
!
!        ...triangular face node
            case(MDLT)
               call nr_face_sons(MDLT,NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  !nods = NODES(nod)%sons(is)
                  nods = Son(nod,is)
                  select case(NODES(nods)%ntype)
                  case(MDLT)
                     nord = max(nord,NODES(nods)%visit)
                  case(MDLQ)
                     call decode(NODES(nods)%visit, nordh,nordv)
                     nord = max(nord,nordv)
                  end select
               enddo
!
!           ...modify the node order
               ! write(*,*) 'call nodmod: nod, nord = ', nod, nord
               call nodmod(nod, nord)
!
!           ...communicate the new order to the (inactive) sons
               call nr_sons(MDLT,NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  !nods = NODES(nod)%sons(is)
                  nods = Son(nod,is)
                  select case(NODES(nods)%ntype)
                  case(MDLT); NODES(nods)%order = nord
                  case(MDLQ); NODES(nods)%order = nord*10+nord
                  end select
               enddo
!
!        ...rectangular face node
            case(MDLQ)
               call decode(nord, nordh,nordv)
               call nr_face_sons(MDLQ,NODES(nod)%ref_kind, nrsons)
               do is=1,nrsons
                  !nods = NODES(nod)%sons(is)
                  nods = Son(nod,is)
                  call decode(NODES(nods)%visit, nordhs,nordvs)
                  nordh = max(nordh,nordhs); nordv = max(nordv,nordvs)
               enddo
!
!           ...modify the node order
               nordhv = max(nordh,nordv)
               nord = nordhv*10+nordhv
               ! nord = nordh*10+nordv
               call nodmod(nod, nord)
!
!           ...communicate the new order to the (inactive) mid-face sons
               do is=1,nrsons
                  !nods = NODES(nod)%sons(is)
                  nods = Son(nod,is)
                  NODES(nods)%order = nord
               enddo
!
!           ...communicate the new order to the (inactive) mid-edge sons
               select case(NODES(nod)%ref_kind)
               case(11)
                  do i=1,4
                     !nods = NODES(nod)%sons(nrsons+i)
                     nods = Son(nod,nrsons+i)
                     select case(i)
                     case(1,3); NODES(nods)%order = nordv
                     case(2,4); NODES(nods)%order = nordh
                     end select
                  enddo
               case(10)
                  !nods = NODES(nod)%sons(nrsons+1)
                  nods = Son(nod,nrsons+1)
                  NODES(nods)%order = nordv
               case(01)
                  !nods = NODES(nod)%sons(nrsons+1)
                  nods = Son(nod,nrsons+1)
                  NODES(nods)%order = nordh
               end select
            end select
!
!     ...unrefined active node
         else
!
!        ...modify the node order
            nord = NODES(nod)%visit
            ! write(*,*) 'call nodmod: nod, nord = ', nod, nord
            call nodmod(nod, nord)
!
!     ...if a refined node
         endif
!
!  ...if an active node
      endif
!
!..end of loop through nodes
   enddo
!
!..reset visitation flags
   call reset_visit
!
end subroutine enforce_max_rule
!
!
!----------------------------------------------------------------------
!> @brief Saves max order in Nod's visit flag
!> @param[in]  Nod   - Node number
!> @param[in]  Nord  - Minimum order requested for the node
!> @date Feb 2023
!----------------------------------------------------------------------
subroutine save_max_order(Nod,Nord)
!
   use data_structure3D
   implicit none
!
   integer, intent(in) :: Nod,Nord
!
   integer :: nordh,nordh1,nordh2
   integer :: nordv,nordv1,nordv2
   integer :: nordhv
!
   if (NODES(Nod)%visit.eq.0) then
      NODES(nod)%visit= Nord
   else
      select case(NODES(Nod)%ntype)
      case(MEDG,MDLT)
         NODES(nod)%visit = max(NODES(nod)%visit,Nord,NODES(nod)%order)
         ! write(*,*) 'save_max_order: Nod, NODES(nod)%visit = ', Nod,NODES(nod)%visit
      case(MDLQ)
         call decode(NODES(nod)%order, nordh,nordv)
         call decode(Nord, nordh1,nordv1)
         call decode(NODES(nod)%visit, nordh2,nordv2)
         nordh = max(nordh,nordh1,nordh2); nordv = max(nordv,nordv1,nordv2)
         ! NODES(nod)%visit = nordh*10+nordv
!     ...enforce isotropic order for the face (TODO: why?)
         nordhv = max(nordh,nordv)
         NODES(nod)%visit = nordhv*10+nordhv
      end select
   endif
!
end subroutine save_max_order
