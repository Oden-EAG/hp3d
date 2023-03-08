      subroutine semrhs (Amdest, Elrhs, Frhs)
c
#include "syscom.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: Assembles the rhs's into the front
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Amdest:  dof destination vectors
c I   Elrhs :  element rhs
c IO  Frhs  :  front assembled rhs
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
#include "surfsc1.blk"
#include "surfsc2.blk"
c
      real(8) :: Amdest(*),Elrhs(*),Frhs(*)
c
      n = 1
c
c loop over the number of rhs's
c
      do 70 in=1,NRHS
c     ---------------
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c loop thru the dof in the element
c
         do 50 i=1,NDOFM
c        ---------------
c pull the dof destination vector
c
            md = Amdest(i)
c
cwb >
c determine the position in the front (function locr)
cwb   (**note: we may wish to inline this code)
c
cwb             j = locr(in,md)
            j = (in-1)*MFW + md
cwb <
c
            Frhs(j) = Frhs(j) + Elrhs(n)
            n = n + 1
   50    continue
c
   70 continue
c
c
      end subroutine semrhs
