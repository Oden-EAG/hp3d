      subroutine zerou (Flhs)
c
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: If necessary, zero out for new unsymmetric lhs equations
c            in the front
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c
c IO  Flhs  : the space which holds the lhs equations in the front
c          x  note: for unsymmetric, lhs is stored as:
c               a11 a12 a13 ...
c               a21 a22 a23 ...  ==>  [a11,a12,a13,..., a21,a22,a23,...]
c               a31 a32 a33 ...
c                :   :   :
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
      use surfsc2
#include "syscomz.blk"
      dimension Flhs(*)
c
      data zero /0.0d0/
c
c
c if the new frontwidth is equal to the old frontwidth
c   then there is nothing to zero
c
      if(NFW .eq. LFW) return
c
      mi = LFW*NFW + 1
      mj = NFW*NFW
      mk = LFW*LFW + 1
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
      do 10 i = mi,mj
        Flhs(i) = zero
   10 continue
c
      if (LFW .eq. 0) return
c
      mj = NFW - LFW
c
      do 50 i = 1,LFW
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
         do 20 j = 1,mj
            mi = mi - 1
            Flhs(mi) = zero
   20    continue
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c also shove to the end the RHS elimination terms
c
         do 30 j = 1,LFW
            mi = mi - 1
            mk = mk - 1
            Flhs(mi) = Flhs(mk)
   30    continue
c
   50 continue
c
c
      end subroutine zerou
