!-------------------------------------------------------------------------------
!> @brief      routine determines BC flags for faces of an element
!!
!! @param[in]  Mdle - middle node
!! @param[out] Ibc  - BC flags; for each variable component ivar supported by
!!                    the element (and its initial mesh ancestor), and iface;
!!             Ibc(iface,ivar) = 0 if the face is interior to the domain,
!!                             = the corresponding flag (1-9) for the element
!!                               ancestor and its face containing the face
!!                               of 'Mdle'.
!> @date       Feb 2023
!-------------------------------------------------------------------------------
!
subroutine find_bc(Mdle, Ibc)
!
   use element_data
   use data_structure3D
   implicit none
!
!..Arguments
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Ibc(6,NRINDEX)
!
!..Locals
   integer :: nodesl(27),norientl(27)
   integer :: ibc_iel(6)
   integer :: nrve,nrf
   integer :: nod,nfath,iel,iface,nrve_iel,nrf_iel,loc,ivar,nvar
!
#if HP3D_DEBUG
   integer :: iprint
   iprint=0
#endif
!
!-------------------------------------------------------------------------------
!
   Ibc = 0
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7010) Mdle
   7010 format(' find_bc: Mdle = ',i8)
      endif
#endif
!
!..determine the initial mesh element ancestor
   nod = Mdle
   do while(nod.gt.0)
     nod = NODES(nod)%father
   enddo
   iel = -nod
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,7020) iel
 7020 format('find_bc: iel = ',i5)
   endif
#endif
!
!..total number of components supported by the initial mesh element
   nvar = ubound(ELEMS(iel)%bcond,1)
!
!..initial mesh element number of vertices + edges, number of faces
   nrve_iel = Nvert(ELEMS(iel)%etype) + Nedge(ELEMS(iel)%etype)
   nrf_iel = Nface(ELEMS(iel)%etype)
!
!..get the element nodes
   call elem_nodes(Mdle, nodesl,norientl)
!
!..number of the element vertices and edges combined
   nrve = Nvert(NODES(Mdle)%ntype) + Nedge(NODES(Mdle)%ntype)
!
!..number of the element faces
   nrf = Nface(NODES(Mdle)%ntype)
!
!..loop through the faces of the element
   do iface=1,nrf
!
!  ...pick up the face node
      nod = nodesl(nrve+iface)
!
!  ...go up the nodal tree
      do
         nfath = NODES(nod)%father
         if (nfath.lt.0) then
!
!        ...initial mesh element
            call locate(nod,ELEMS(iel)%nodes(nrve_iel+1:nrve_iel+nrf_iel),nrf_iel, loc)
            do ivar=1,nvar
               call decodg(ELEMS(iel)%bcond(ivar),10,nrf_iel, ibc_iel)
#if HP3D_DEBUG
               if (iprint.eq.1) then
             7030 format('find_bc: BC FOR iel=',i8,' AND ivar=',i2,': ',6i2)
                  write(*,7030) ivar,ibc_iel(1:nrf_iel)
               endif
#endif
               Ibc(iface,ivar) = ibc_iel(loc)
            enddo
            goto 10
         else
            select case(NODES(nfath)%ntype)
               case(MDLT,MDLQ)
!
!           ...a mid-face node, continue up the tree
               nod = nfath
               case default
!
!           ...a middle node, quit, the face is in the interior of the domain
               goto 10
            end select
         endif
!  ...end loop through nodal tree
      enddo
   10 continue
!..end loop through element faces
   enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
         do ivar=1,nvar
       7100 format('          ivar=',i2,',  Ibc = ',6(i1,2x))
            write(*,7100) ivar,Ibc(1:nrf,ivar)
         enddo
      endif
#endif
!
!
end subroutine find_bc
