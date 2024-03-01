!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! function to determine the position in the current front for
!  a symmetric lhs array element k(i,j)
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
      function locls (Iin,Jin)
      implicit none
!
      integer :: locls
      integer :: Iin,Jin
      integer :: ii,jj,k,l
!
      ii = abs(Iin)
      jj = abs(Jin)
!
      k = max0(ii,jj)
      l = min0(ii,jj)
!
      locls = k*(k-1)/2 + l
!
!
      end function locls
