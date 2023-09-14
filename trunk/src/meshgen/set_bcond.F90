!-------------------------------------------------------------------------------
!> @brief      Set boundary condition flags for a physics attribute component
!!             of all initial mesh elements if the face domain number matches
!!             the requested boundary domain ID.
!!
!> @param[in]  Dom  - boundary domain ID where to set BC on
!> @param[in]  Attr - physics attribute to set BC for: 1,...,NR_PHYSA
!> @param[in]  Comp - attribute component to set BC for: 1,...,NR_COMP(Attr)
!> @param[in]  Flag - BC flag: 1,...,9
!!
!> @date       Sep 2023
!-------------------------------------------------------------------------------
subroutine set_bcond(Dom,Attr,Comp,Flag)
!
   use data_structure3D, only: NRELIS
!
   implicit none
!
   integer, intent(in) :: Dom,Attr,Comp,Flag
!
   integer :: iel
!
!$OMP PARALLEL DO
!..loop through initial mesh elements
   do iel=1,NRELIS
      call set_bcond_elem(iel,Dom,Attr,Comp,Flag)
   enddo
!$OMP END PARALLEL DO
!
end subroutine set_bcond


!-------------------------------------------------------------------------------
!> @brief      Set boundary condition flags for a physics attribute component
!!             of an initial mesh element if the face domain number matches
!!             the requested boundary domain ID.
!!
!> @param[in]  Iel  - initial mesh element: 1,...,NRELIS
!> @param[in]  Dom  - boundary domain ID where to set BC on
!> @param[in]  Attr - physics attribute to set BC for: 1,...,NR_PHYSA
!> @param[in]  Comp - attribute component to set BC for: 1,...,NR_COMP(Attr)
!> @param[in]  Flag - BC flag: 1,...,9
!!
!> @date       Sep 2023
!-------------------------------------------------------------------------------
subroutine set_bcond_elem(Iel,Dom,Attr,Comp,Flag)
!
   use data_structure3D, only: ELEMS,NRELIS
   use GMP             , only: TRIANGLES,RECTANGLES,HEXAS,TETRAS,PRISMS,PYRAMIDS
   use node_types
   use physics
!
   implicit none
!
   integer, intent(in) :: Iel,Dom,Attr,Comp,Flag
!
!..auxiliary array with BC flags for each face
   integer :: ibc(1:6)
!
   integer :: etype,ifc,idom,index,neig,nelem,nrect,ntria
!
!..verifying input arguments
   if (Iel.lt.1 .or. Iel.gt.NRELIS) then
      write(*,1000) 'Iel',Iel
      stop
   endif
   if (Attr.lt.1 .or. Attr.gt.NR_PHYSA) then
      write(*,1000) 'Attr',Attr
      stop
   endif
   if (Comp.lt.1 .or. Comp.gt.NR_COMP(Attr)) then
      write(*,1000) 'Comp',Comp
      stop
   endif
   if (Flag.lt.1 .or. Flag.gt.9) then
      write(*,1000) 'Flag',Flag
      stop
   endif
   1000 format('set_bcond_elem: invalid input: ',A,' = ',I9)
!
!..set BC flags on all exterior faces of this element
!  with the correct boundary domain ID
   ibc(1:6) = 0
!
   etype = ELEMS(Iel)%etype
   do ifc=1,NFACE(etype)
      neig = ELEMS(Iel)%neig(ifc)
!  ...exterior face (no neighbor)
      if (neig .eq. 0) then
!     ...get element number in GMP data structure
         nelem = ELEMS(Iel)%GMPblock / 10
!     ...get domain number of this face
         select case(etype)
            case(BRIC)
               nrect = HEXAS(nelem)%FigNo(ifc) / 10
               idom = RECTANGLES(nrect)%Domain
            case(TETR)
               ntria = TETRAS(nelem)%FigNo(ifc) / 10
               idom = TRIANGLES(ntria)%Domain
            case(PRIS)
               select case(ifc)
                  case(1,2)
                     ntria = PRISMS(nelem)%FigNo(ifc) / 10
                     idom = TRIANGLES(ntria)%Domain
                  case(3,4,5)
                     nrect = PRISMS(nelem)%FigNo(ifc) / 10
                     idom = RECTANGLES(nrect)%Domain
               end select
            case(PYRA)
               select case(ifc)
                  case(1)
                     nrect = PYRAMIDS(nelem)%FigNo(ifc) / 10
                     idom = RECTANGLES(nrect)%Domain
                  case(2,3,4,5)
                     ntria = PYRAMIDS(nelem)%FigNo(ifc) / 10
                     idom = TRIANGLES(ntria)%Domain
               end select
         end select
!     ...set BC flag if the face domain ID matches the input boundary domain ID
         if (idom .eq. Dom) ibc(ifc) = Flag
      endif
!..end loop over element faces
   enddo
!
!..Check if BC flag array has been allocated
   if (.not. associated(ELEMS(Iel)%bcond)) then
      allocate(ELEMS(Iel)%bcond(NRINDEX))
   endif
!
!..Encode face BC into a single BC flag
   call attr_to_index(Attr,Comp, index)
   call encodg(ibc(1:6),10,6, ELEMS(Iel)%bcond(index))
!
end subroutine set_bcond_elem
