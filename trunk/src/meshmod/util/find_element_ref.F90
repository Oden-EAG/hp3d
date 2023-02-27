!> @brief      modify the refinement to accommodate existing refinements
!!
!> @param[in]  Ntype - middle node type
!> @param[in]  Kref  - refinement kind of middle node
!> @param[in]  Kreff - existing face refinements
!!
!> @param[out] Krefm - modified (upgraded) refinement
!!
!> @date       Feb 2023
subroutine find_element_ref(Ntype,Kref,Kreff, Krefm)
  use element_data
  use refinements
  implicit none
  ! ** Arguments
  integer,               intent(in)  :: Ntype
  integer,               intent(in)  :: Kref
  integer, dimension(6), intent(in)  :: Kreff
  integer,               intent(out) :: Krefm
  ! ** Locals
  integer, dimension(6) :: kreff_trial
  integer :: ipass, isum, iface, i, j, kref_trial
  !
#if DEBUG_MODE
   integer :: iprint = 0
#endif
  !
  !-------------------------------------------------------
  !
#if DEBUG_MODE
  if (iprint.ge.1) then
      write(*,7000) S_Type(NType),Kref,Kreff(1:6)
  endif
#endif
  !
  7000 format(' find_element_ref: Type,Kref,Kreff = ',a4,',',i3,',',6(2x,i2))
  7001 format(' find_element_ref: Type,Kref,Krefm,Kreff = ',   &
                                  a4,',',i3,',',i3,',',6(2x,i2))
  !
  Krefm = -1
  ! be careful, the search direction is from aniso refinement
  ! to minimize the unwanted refinements
  do i=nr_ref(Ntype),1,-1
     kref_trial = kref_kind(i,Ntype)
     call check_ref(Ntype,Kref,kref_trial, ipass)
     !
     ! if not pass, try rest
     if (ipass.eq.0) cycle
     !
     ! if pass element level refinement, check face
     call find_face_ref_flags(Ntype,kref_trial, kreff_trial)
     !
     isum = 0
     do iface=1,nface(Ntype)
        j = nvert(Ntype) + nedge(Ntype) + iface
        call check_ref(Type_nod(Ntype,j), Kreff(iface), &
                       kreff_trial(iface), ipass)
        isum = isum + ipass
     enddo
     !
     ! pass all test, loop out with modified one
     if (isum.eq.nface(Ntype)) then
        Krefm = kref_trial
        exit
     endif
  enddo
  !
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7001) S_Type(Ntype),Kref,Krefm,Kreff(1:nface(Ntype))
     call pause
  endif
#endif
  !
end subroutine find_element_ref
