!> Purpose : find out the minimum refinement for closing
!! @param[in]  Type  - middle node type
!! @param[in]  Kreff - refinement flags for faces
!! @param[in]  Krefe - refinement flags for edges
!! @param[out] Kref  - refinement flag for the middle node
subroutine find_element_closing_ref(Type,Kreff,Krefe, Kref)
  use refinements
  use element_data
  implicit none
  ! ** Arguments
  !-------------------------------------------------------
  character(len=4),       intent(in)  :: Type
7001 format('find_element_closing_ref: Type = ',a5)
  integer, dimension(6),  intent(in)  :: Kreff
  integer, dimension(12), intent(in)  :: Krefe
  integer,                intent(out) :: Kref
  ! ** Locals
  !-------------------------------------------------------
  integer, dimension(6)  :: Kreff_trial
  integer, dimension(12) :: Krefe_trial
  
  integer :: iprint, i, j, isum, ipass, iref, kref_trial
  !-------------------------------------------------------
  iprint=0

  ! check anisotropic refinement first
  do iref=nr_ref(Type),1,-1
     kref_trial = kref_kind(iref, Type)

     ! Step 0 : trial refinement
     !~~~~~~~~~~~~~~~~~~~~~~~~~~
     call find_face_ref_flags(Type,kref_trial, kreff_trial)
     call find_edge_ref_flags(Type,kref_trial, krefe_trial)

     ! test
     isum = 0

     ! Step 1 : check face
     !~~~~~~~~~~~~~~~~~~~~~
     do i=1,nface(Type)
        j = nvert(Type) + nedge(Type) + i
        call check_ref(Type_nod(Type, j), &
                            Kreff(i),kreff_trial(i), ipass)
        isum = isum + ipass
     enddo

     ! Step 2 : check edge
     !~~~~~~~~~~~~~~~~~~~~
     do i=1,nedge(Type)
        if (krefe_trial(i).lt.Krefe(i)) then
           ipass = 0
        else
           ipass = 1
        endif
        isum = isum + ipass
     enddo

     ! Step 3 : if all pass, set the kref
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if (isum.eq.(nface(Type)+nedge(Type))) then
        Kref = kref_trial; exit
     endif
  enddo

  if (iprint.eq.1) then
     write(*,7001) Type
     write(*,7002) Kreff(1:nface(Type))
7002 format('                          Face flags = ',6i3)
     write(*,7003) Krefe(1:nedge(Type))
7003 format('                          Edge flags = ',12i3)
     write(*,7004) Kref
7004 format('                          Element flag = ',i3)
  endif

end subroutine find_element_closing_ref
