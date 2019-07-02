!> Purpose : modify the refinement to accommodate existing refinements
!! @param[in]  Type  - middle node type
!! @param[in]  Kref  - refinement kind of middle node
!! @param[in]  Kreff - existing face refimenets
!!
!! @param[out] Krefm - modified refinement
subroutine find_element_ref(Type,Kref,Kreff, Krefm)
  use element_data
  use refinements
  implicit none
  ! ** Arguments
#ifdef _PYHP3D
  integer,      intent(in)  :: Type
 7000 format(' find_element_ref: Type,Kref,Kreff = ',i2,4x,i2,4x,6(i2,2x))
7001 format('find_element_ref: Type,Kref,Krefm,Kreff = ', &
          i2,2x,i2,2x,i2,6x,6(2x,i2))
#else
  character(len=4),      intent(in)  :: Type
 7000 format(' find_element_ref: Type,Kref,Kreff = ',a4,4x,i2,4x,6(i2,2x))
7001 format('find_element_ref: Type,Kref,Krefm,Kreff = ', &
          a5,2x,i2,2x,i2,6x,6(2x,i2))
#endif
  integer,               intent(in)  :: Kref
  integer, dimension(6), intent(in)  :: Kreff
  integer,               intent(out) :: Krefm
  ! ** Locals
  integer, dimension(6) :: kreff_trial
  integer :: iprint, ipass, isum, iface, i, j, kref_trial
  !-------------------------------------------------------

  iprint=0

  if (iprint.ge.1) then
      write(*,7000)Type,Kref,Kreff(1:6)
  endif

  Krefm = -1
  ! be careful, the search direction is from aniso refinement
  ! to minimize the unwanted refinements
  do i=nr_ref(Type),1,-1
     kref_trial = kref_kind(i, Type)
     call check_ref(Type,Kref,kref_trial, ipass)
     
     ! if not pass, try rest
     if (ipass.eq.0) cycle

     ! if pass element level refinement, check face
     call find_face_ref_flags(Type,kref_trial, kreff_trial)

     isum = 0
     do iface=1,nface(Type)
        j = nvert(Type) + nedge(Type) + iface
        call check_ref(Type_nod(Type, j), Kreff(iface), &
                       kreff_trial(iface), ipass)
        isum = isum + ipass
     enddo

     ! pass all test, loop out with modified one
     if (isum.eq.nface(Type)) then
        Krefm = kref_trial
        exit
     endif
  enddo
!
  if (iprint.eq.1) then
     write(*,7001) Type,Kref,Krefm,Kreff(1:nface(Type))
     call pause
  endif

end subroutine find_element_ref
