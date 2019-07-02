      subroutine elmrhs (Aid, Inc, Frhs, Ubuf, Bbuf)
c
#include "syscom.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c Elimination of rhs's for equation (id)
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c
c I   Aidx   transfer loc for active dof in the front after eliminating
c             this is the equation to be eliminated
c
c I   Inc    for unsymmetric this is the Increment to shift for picking
c             the correct entries in Ubuf()
c
c IO  Frhs   the assembled front rhs
c
c IO  Ubuf   the buffer to hold eliminated lhs equations
c IO  Bbuf   the buffer to hold eliminated rhs equations
c
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
      dimension Frhs(*), Ubuf(*), Bbuf(*)
cwb >
      data sml2 /1.e-15/
cwb <

      id = Aid
      idm = id - 1
      idp = id + 1
      im = 0
c
c perform the rhs elimination:
c-----------------------------
c   f(i)' = f(i) - [k(i,m)/k(m,m)]*f(m)
c
c  where m is the row being eliminated
c
c
c loop over the number of rhs in the model
c
      do 60 in = 1,NRHS
c     -----------------
c
         iuu = 1
c
c pull f(m)
c
         s = Frhs(im+id)
c
c ***note: we load the rhs elimination buffer with the last  f(i)'
c           (ie: f(i)" = (f(i)'-[k(i,m)/k(m,m)]*f(m))....)
c
         Bbuf(in) = s
c        ----------
cwb >
c dont operate on zero elements
c -----------------------------
         if (dabs(s) .le. sml2) then
            iuu = iuu + idm*Inc
            go to 25
         endif
cwb <
c
         ii1 = im + 1
         ii2 = im + idm
c
c loop over the front dof, up to the dof we are eliminating
c note: u() contains already the entries [k(i,m)/k(m,m)]
c and eliminate
c     ---------
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c
         do 20 ii = ii1,ii2
            Frhs(ii) = Frhs(ii) - s*Ubuf(iuu)
            iuu = iuu + Inc
   20    continue
c
   25    if (ISYM.eq.1 .or. ISYM.eq.4) iuu = iuu + 1
c
cwb          if (idp .gt. NFW) go to 50
c
c loop thru the dof below the elimination dof in the front
c----------------------------------------------------------
c **note: we are suffling the equations forward as we perform the calcs
c
        if (idp .le. NFW) then
c
            ii1 = im + idp
            ii2 = im + NFW
c
cwb >
c dont operate on zero elements
c -----------------------------
            if (dabs(s) .le. sml2) then
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c
               do 27 ii = ii1,ii2
                  Frhs(ii-1) = Frhs(ii)
   27          continue
c
               iuu = iuu + (NFW - idp + 1)*Inc
            else
cwb <
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c
               do 30 ii = ii1,ii2
                  Frhs(ii-1) = Frhs(ii) - s*Ubuf(iuu)
                  iuu = iuu + Inc
   30          continue
            endif
         endif
c
         im = im + MFW
   60 continue
c
      return
      end
