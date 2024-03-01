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
      real(8) :: Flhs(*)
!
      real(8), parameter :: dzero = 0.0d0
!
      integer :: i,mi,mj
!
!  ...if the new frontwidth is equal to the old frontwidth
!     then there is nothing to zero
!
      if (NFW .eq. LFW) return
!
! set the zeroing bounds
!
      mi = (LFW*(LFW+1))/2 + 1
      mj = (NFW*(NFW+1))/2
!
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!
      do 10 i = mi,mj
         Flhs(i) = dzero
   10 continue
!
!
      end subroutine zeros
