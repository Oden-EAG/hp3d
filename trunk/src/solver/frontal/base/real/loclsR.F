      function locls (Iin,Jin)
c
#include "syscom.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c function to determine the position in the current front for
c  a symmetric lhs array element k(i,j)
c    ie: the contribution for dof destvec i by dof destvec j
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Iin : row number into the current front for k(i,j)
c I   Jin : column number into the current front for k(i,j)
c*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
c AUTHOR:                             LATEST REVISION:
c++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
c NAMING CONVENTIONS:
c     AAAAAAAA    Variables in COMMON & PARAMETERS
c     Aaaaaaaa    Variables as ARGUMENTS
c     aaaaaaaa    LOCAL Variables
c         7xxx    FORMAT Statements
c         9xxx    ERROR Handling
c+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
c
c
      ii = abs(Iin)
      jj = abs(Jin)
c
      k = max0(ii,jj)
      l = min0(ii,jj)
c
      locls = k*(k-1)/2 + l
c
c
      end function locls
