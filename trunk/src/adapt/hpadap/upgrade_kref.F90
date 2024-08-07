!-----------------------------------------------------------------------
!
!   routine name      - upgrade_kref
!
!-----------------------------------------------------------------------
!
!   latest revision    - June 2023
!
!   purpose            - upgrades the refinement flag to
!                        provide a union of the intended flag and the existing flag
!
!   arguments :
!     in:        kref: existing flag
!     out:       krefm (in-out)
!                in:  intended flag
!                out: union
!
!-----------------------------------------------------------------------
subroutine upgrade_kref(kref, krefm)
   implicit none
   integer, intent(in)    :: kref
   integer, intent(out) :: krefm
!
   integer :: i,ii(3),jj(3)
!
   call ddecode(kref,  ii(1),ii(2),ii(3))
   call ddecode(krefm, jj(1),jj(2),jj(3))
!
   do i=1,3
      jj(i) = min(1, ii(i)+jj(i))
   enddo
!
   krefm = jj(1)*100 + jj(2)*10 + jj(3)
!  
end subroutine upgrade_kref