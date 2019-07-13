!-------------------------------------------
! routine emulates obsolete Fortran pause
!-------------------------------------------
subroutine pause
   write(*,*) 'PAUSE - type anything'
   read(*,*)
   return
end subroutine
