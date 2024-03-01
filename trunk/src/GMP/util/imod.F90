!----------------------------------------------------------------------
!> @brief function computes the modulus of i/m . Differently than
!!           the built-in function "mod", for i > 0, function imod sets
!!           representative for equivalence class of 0 equal to m,
!!           instead of 0; example :
!!             imod(3,3) = 3     ;     mod(3,3) = 0
!!
!> @date Nov 12
!----------------------------------------------------------------------
integer function imod(i,m)
      implicit none
      integer,intent(in) :: i,m
      imod=i-(i-1)/m*m
endfunction imod
!
!
! For backward compatibility:
integer function mod3(i)
      implicit none
      integer,intent(in) :: i
      integer,external :: imod
      mod3=imod(i,3)
endfunction mod3
!
!
integer function mod4(i)
      implicit none
      integer,intent(in) :: i
      integer,external :: imod
      mod4=imod(i,4)
endfunction mod4
!
!
integer function my_mod(i,m)
      implicit none
      integer,intent(in) :: i,m
      integer,external :: imod
      my_mod=imod(i,m)
endfunction my_mod
