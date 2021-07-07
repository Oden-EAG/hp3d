!> Purpose : modify the refinement flag for an element to accommodate existing face refinements
!! @param[in]  Type  - middle node type
!! @param[in]  Kref  - refinement kind of middle node
!! @param[in]  Kreff - existing face refinements
!!
!! @param[out] Krefm - modified (upgraded) refinement
subroutine find_element_ref(Type,Kref,Kreff, Krefm)
  use element_data
  use refinements
  implicit none
  ! ** Arguments
  character(len=4),      intent(in)  :: Type
  integer,               intent(in)  :: Kref
  integer, dimension(6), intent(in)  :: Kreff
  integer,               intent(out) :: Krefm
  ! ** Locals
  integer, dimension(6) :: kreff_trial
  integer :: iprint, ipass, isum, iface, i, j, kref_trial
  !-------------------------------------------------------
  !
  iprint=1
#if DEBUG_MODE
  if (iprint.ge.1) then
      write(*,7000) Type,Kref,Kreff(1:6)
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
  do i=nr_ref(Type),1,-1
     kref_trial = kref_kind(i, Type)
     call check_ref(Type,Kref,kref_trial, ipass)
     !
     ! if has not passed, exit
     if (ipass.eq.0) cycle
     !
     ! if it has passed element refinement test, check faces
     call find_face_ref_flags(Type,kref_trial, kreff_trial)
     !
     isum = 0
     do iface=1,nface(Type)
        j = nvert(Type) + nedge(Type) + iface
        call check_ref(Type_nod(Type, j), Kreff(iface), &
                       kreff_trial(iface), ipass)
        isum = isum + ipass
     enddo
     !
     ! has passed the whole test, loop out with the modified refinement flag
     if (isum.eq.nface(Type)) then
        Krefm = kref_trial
        exit
     endif
  enddo
  !
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7001) Type,Kref,Krefm,Kreff(1:nface(Type))
     call pause
  endif
#endif
  !
end subroutine find_element_ref
