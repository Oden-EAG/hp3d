c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c function to determine the position in the current front for
c  a unsymmetric lhs array element k(i,j)
c    ie: the contribution for dof destvec i by dof destvec j
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Iin : row number into the current front for k(i,j)
c I   Jin : column number into the current front for k(i,j)
c*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
c LATEST REVISION: Feb 2023
c++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
c NAMING CONVENTIONS:
c     AAAAAAAA    Variables in COMMON & PARAMETERS
c     Aaaaaaaa    Variables as ARGUMENTS
c     aaaaaaaa    LOCAL Variables
c         7xxx    FORMAT Statements
c         9xxx    ERROR Handling
c+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
c
      function loclu (Iin,Jin)
c
#include "surfsc2.blk"
c
      integer :: loclu
      integer :: Iin,Jin
      integer :: ii,jj
c
      ii = abs(Iin)
      jj = abs(Jin)
      loclu = jj + (ii-1)*NFW
c
      end function loclu
