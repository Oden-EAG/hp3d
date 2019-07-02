      subroutine zeros (Flhs)
c
#include "syscomz.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: If necessary, zero out for new symmetric lhs equations
c            in the front
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c
c IO  Flhs  : the space which holds the lhs equations in the front
c             note: for symmetric, lhs is stored as:  *** by column ***
c               a11 a12 a13 ...
c                   a22 a23 ...  ==>   [a11, a12,a22, a13,a23,a33, ...]
c                       a33 ...
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
#include "surfsc2.blk"
c
      dimension Flhs(*)
c
      data zero /0.0d0/
c
c if the new frontwidth is equal to the old frontwidth
c   then there is nothing to zero
c
      if (NFW .eq. LFW) return
c
c set the zeroing bounds
c
      mi = (LFW*(LFW+1))/2 + 1
      mj = (NFW*(NFW+1))/2
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
      do i = mi,mj
        Flhs(i) = zero
      enddo
c
      return
      end
