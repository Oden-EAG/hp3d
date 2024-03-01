!-------------------------------------------------------------------------------
!> @brief       Set physical attributes to be supported by the all elements
!> @param[in]   NrAttr   - number of supported attributes: 1,...,NR_PHYSA
!> @param[in]   AttrList - list of supported attributes:
!!                        AttrList(i) \in {1,...,NR_PHYSA}, i = 1,...,NrAttr
!> @date        Sep 2023
!-------------------------------------------------------------------------------
subroutine set_attr(NrAttr,AttrList)
!
   use data_structure3D, only: NRELIS
!
   implicit none
!
   integer, intent(in) :: NrAttr
   integer, intent(in) :: AttrList(NrAttr)
!
   integer :: iel
!
!$OMP PARALLEL DO
!..loop through initial mesh elements
   do iel=1,NRELIS
      call set_attr_elem(iel,NrAttr,AttrList)
   enddo
!$OMP END PARALLEL DO
!
end subroutine set_attr


!-------------------------------------------------------------------------------
!> @brief       Set physical attributes to be supported by the element
!> @param[in]   Iel      - initial mesh element: 1,...,NRELIS
!> @param[in]   NrAttr   - number of supported attributes: 1,...,NR_PHYSA
!> @param[in]   AttrList - list of supported attributes:
!!                        AttrList(i) \in {1,...,NR_PHYSA}, i = 1,...,NrAttr
!> @date        Sep 2023
!-------------------------------------------------------------------------------
subroutine set_attr_elem(Iel,NrAttr,AttrList)
!
   use data_structure3D, only: ELEMS,NRELIS
   use physics         , only: NR_PHYSA,PHYSA
!
   implicit none
!
   integer, intent(in) :: Iel,NrAttr
   integer, intent(in) :: AttrList(NrAttr)
!
   integer :: iattr,attr
!
!..verifying input arguments
   if (Iel.lt.1 .or. Iel.gt.NRELIS) then
      write(*,1000) 'Iel',Iel
      stop
   endif
   if (NrAttr.gt.NR_PHYSA) then
      write(*,1000) 'NrAttr',NrAttr
      stop
   endif
   1000 format('set_attr: invalid input: ',A,' = ',I9)
!
!..set physical attributes supported by the element
   ELEMS(Iel)%nrphysics=NrAttr
   allocate(ELEMS(Iel)%physics(NrAttr))
!
   do iattr = 1,NrAttr
      attr = AttrList(iattr)
      if (attr.lt.1 .or. attr.gt.NR_PHYSA) then
         write(*,1000) 'attr',attr
         stop
      endif
      ELEMS(Iel)%physics(iattr) = PHYSA(attr)
   enddo
!
end subroutine set_attr_elem
