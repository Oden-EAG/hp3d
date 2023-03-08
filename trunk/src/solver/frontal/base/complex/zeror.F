       subroutine zeror (Frhs)
c
#include "syscomz.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: If necessary, zero out for new rhs equations in the front
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c IO  Frhs  : the space which holds the rhs equations in the front
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
#include "surfsc1.blk"
#include "surfsc2.blk"
c
      dimension Frhs(*)
c
      data zero /0.0d0/
c
      iprint=0
c
c if the new frontwidth is equal to the old frontwidth
c   then there is nothing to zero
c
      if (LFW .eq. NFW) return
c
      if (iprint.eq.1) then
        write(nfsout,*) 'ZEROR: NRHS,MFW,NFW,LFW = ',NRHS,MFW,NFW,LFW
      endif
      do 120 in  =1,NRHS
c
         ia = (in-1) * MFW
         m = LFW + 1
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
         do 20 i = m,NFW
            Frhs(ia+i) = zero
20       continue
c
120   continue
c
c
      end subroutine zeror
