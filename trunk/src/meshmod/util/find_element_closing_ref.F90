!> @brief      find out the minimum refinement for closing
!> @param[in]  Ntype - middle node type
!> @param[in]  Kreff - refinement flags for faces
!> @param[in]  Krefe - refinement flags for edges
!> @param[out] Kref  - refinement flag for the middle node
!> @date       Feb 2023
subroutine find_element_closing_ref(Ntype,Kreff,Krefe, Kref)
  use refinements
  use element_data
  implicit none
  ! ** Arguments
  !-------------------------------------------------------
  integer,                intent(in)  :: Ntype
7001 format('find_element_closing_ref: Type = ',a5)
  integer, dimension(6),  intent(in)  :: Kreff
  integer, dimension(12), intent(in)  :: Krefe
  integer,                intent(out) :: Kref
  ! ** Locals
  !-------------------------------------------------------
  integer, dimension(6)  :: Kreff_trial
  integer, dimension(12) :: Krefe_trial

  integer :: i, j, isum, ipass, iref, kref_trial
  !
#if HP3D_DEBUG
  integer :: iprint
  iprint=0
#endif
  !
  !-------------------------------------------------------
  !
  ! check anisotropic refinement first
  do iref=nr_ref(Ntype),1,-1
     kref_trial = kref_kind(iref,Ntype)

     ! Step 0 : trial refinement
     !~~~~~~~~~~~~~~~~~~~~~~~~~~
     call find_face_ref_flags(Ntype,kref_trial, kreff_trial)
     call find_edge_ref_flags(Ntype,kref_trial, krefe_trial)

     ! test
     isum = 0

     ! Step 1 : check face
     !~~~~~~~~~~~~~~~~~~~~~
     do i=1,nface(Ntype)
        j = nvert(Ntype) + nedge(Ntype) + i
        call check_ref(TYPE_NOD(j,Ntype), &
                       Kreff(i),kreff_trial(i), ipass)
        isum = isum + ipass
     enddo

     ! Step 2 : check edge
     !~~~~~~~~~~~~~~~~~~~~
     do i=1,nedge(Ntype)
        if (krefe_trial(i).lt.Krefe(i)) then
           ipass = 0
        else
           ipass = 1
        endif
        isum = isum + ipass
     enddo

     ! Step 3 : if all pass, set the kref
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if (isum.eq.(nface(Ntype)+nedge(Ntype))) then
        Kref = kref_trial; exit
     endif
  enddo
  !
#if HP3D_DEBUG
  if (iprint.eq.1) then
     write(*,7001) S_Type(Ntype)
     write(*,7002) Kreff(1:nface(Ntype))
7002 format('                          Face flags = ',6i3)
     write(*,7003) Krefe(1:nedge(Ntype))
7003 format('                          Edge flags = ',12i3)
     write(*,7004) Kref
7004 format('                          Element flag = ',i3)
  endif
#endif
  !
end subroutine find_element_closing_ref
