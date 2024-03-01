!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! function to determine the position in the current front for
!  a rhs array element f(i)  for rhs number n
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   Nrh  : the rhs number
! I   Iin  : the current equation being assembled into the front
!             (ie: the current dof destination vector)
!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
! LATEST REVISION: Feb 2023
!++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
! NAMING CONVENTIONS:
!     AAAAAAAA    Variables in COMMON & PARAMETERS
!     Aaaaaaaa    Variables as ARGUMENTS
!     aaaaaaaa    LOCAL Variables
!         7xxx    FORMAT Statements
!         9xxx    ERROR Handling
!+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
!
      function locr (Nrh,Iin)
!
      use surfsc2
      implicit none
!
      integer :: locr
      integer :: Nrh,Iin
      integer :: ii
!
      ii = abs(Iin)
      locr = (Nrh-1)*MFW + ii
!
      end function locr
