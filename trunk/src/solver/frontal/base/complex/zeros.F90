!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: If necessary, zero out for new symmetric lhs equations
!            in the front
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
!
! IO  Flhs  : the space which holds the lhs equations in the front
!             note: for symmetric, lhs is stored as:  *** by column ***
!               a11 a12 a13 ...
!                   a22 a23 ...  ==>   [a11, a12,a22, a13,a23,a33, ...]
!                       a33 ...
!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
! LATEST REVISION: Mar 2023
!++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
! NAMING CONVENTIONS:
!     AAAAAAAA    Variables in COMMON & PARAMETERS
!     Aaaaaaaa    Variables as ARGUMENTS
!     aaaaaaaa    LOCAL Variables
!         7xxx    FORMAT Statements
!         9xxx    ERROR Handling
!+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
!
   subroutine zeros (Flhs)
!
      use surfsc2
!
      implicit none
!
      complex(8) :: Flhs(*)
!
      integer :: i,mi,mj
!
      complex(8), parameter :: zero = (0.d0, 0.d0)
!
! if the new frontwidth is equal to the old frontwidth
!   then there is nothing to zero
!
      if (NFW .eq. LFW) return
!
! set the zeroing bounds
!
      mi = (LFW*(LFW+1))/2 + 1
      mj = (NFW*(NFW+1))/2
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
      do i = mi,mj
        Flhs(i) = zero
      enddo
!
!
   end subroutine zeros
