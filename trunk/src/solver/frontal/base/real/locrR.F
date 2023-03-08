c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c function to determine the position in the current front for
c  a rhs array element f(i)  for rhs number n
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Nrh  : the rhs number
c I   Iin  : the current equation being assembled into the front
c             (ie: the current dof destination vector)
c*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
c LATEST REVISION: Mar 2023
c++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
c NAMING CONVENTIONS:
c     AAAAAAAA    Variables in COMMON & PARAMETERS
c     Aaaaaaaa    Variables as ARGUMENTS
c     aaaaaaaa    LOCAL Variables
c         7xxx    FORMAT Statements
c         9xxx    ERROR Handling
c+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
c
      function locr (Nrh,Iin)
c
#include "surfsc2.blk"
c
      integer :: locr
      integer :: Nrh,Iin
      integer :: ii
c
      ii = abs(Iin)
      locr = (Nrh-1)*MFW + ii
c
c
      end function locr
