!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! function to determine the position in the current front for
!  a unsymmetric lhs array element k(i,j)
!    ie: the contribution for dof destvec i by dof destvec j
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   Iin : row number into the current front for k(i,j)
! I   Jin : column number into the current front for k(i,j)
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
      function loclu (Iin,Jin)
!
      use surfsc2
!
      integer :: loclu
      integer :: Iin,Jin
      integer :: ii,jj
!
      ii = abs(Iin)
      jj = abs(Jin)
      loclu = jj + (ii-1)*NFW
!
      end function loclu
