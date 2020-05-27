      subroutine elmsol (Aidx, Ubuf, Bbuf, Xfrnt)
c
#include "syscomz.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c Calculates the solution for one dof specified by (id)
c
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Aidx  : position of the dof being solved for within the current fr
c             (more specifically the relative position in u()
c
c I   Ubuf  : the eliminated lhs buffered equations
c              positioned relative to the this equation
c
c I   Bbuf  : the eliminated rhs buffered equations
c              positioned relative to the this equation
c
c IO  Xfrnt : the solution front
cc*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
c AUTHOR:                             LATEST REVISION:
c++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
c NAMING CONVENTIONS:
c     AAAAAAAA    Variables in COMMON & PARAMETERS
c     Aaaaaaaa    Variables as ARGUMENTS
c     aaaaaaaa    LOCAL Variables
c         7xxx    FORMAT Statements
c         9xxx    ERROR Handling
c+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
c Additional notes:
c for a given element the back-substituion will look something like:
c  (note: this implies subsequent call to sub.elmsol
c
c |k11  -  k13 k14  - |       x3 =  f3/k33
c |k21 k22 k23 k24 k25|  ===> x4 = (f4 - k34*x3)/k44
c | -   -  k33  -   - |       x1 = (f1 - k13*x3 - k14*x4)/k11
c | -   -  k34 k44  - |       x5 = (f5 - k53*x3 - k54*x4 - k51*x1)/k55
c |k51  -  k53 k54 k55|       x2 = (f2 - k23*x3 - k24*x4 - k21*x1 - k25*
c
c where u would be stored as:
c  Ubuf = {..., k21 k22 k23 k24 k25, k51 k53 k54 k55, k11 k13 k14, k34 k
c  Bbuf = {..., f2, f5, f1, f4, f3 }
c  Xfrnt= {..., x2, x5, x1, x4, x3 }
c
c note: that it doesnt look like this forever,( ie: this is really a fro
c        some equations remain till the next element, some are 'forgotte
c        they have already been passed to permanent storage, and the 'fo
c        ones become new destinations in Xfrnt() based on Aidx
c
#include "surfsc1.blk"
#include "surfsc2.blk"
c
      dimension Ubuf(*), Bbuf(*), Xfrnt(*)


c
      data one /1.0/




c
      iprint=0
c
      id = Aidx
      idm = id - 1
      idp = id + 1
      if (iprint.eq.1) then
        write(nfsout,*) 'ELMSOL: Aidx = ',Aidx
      endif
c
c because Ubuf contains slightly different info for symmetric versus uns
c  we must take a different path herein
c
c  for symmetric: Ubuf() contains:
c       [k(1,m)/k(m,m), k(2,m)/k(m,m),..., k(m,m), ...,k(nfw,m)/k(m,m)]
c  for unsymmetric: Ubuf() contains:
c       [k(1,m), k(2,m),...k(nfw,m)]
c --------------------------------------------------------------------
c pull the k(m,m) entry from the Ubuf buffer
c-----------------------------------------
      if (ISYM.eq.1 .or. ISYM.eq.4) then
         f1 = Ubuf(id)
         f2 = one
      else
         f1 = one
         f2 = Ubuf(id)
      endif
c
c loop thru and back-substitute for each rhs
c --------------------------------------------
      do 40 in = 1,NRHS
c     -----------------
         iuu = NFW


c
c ja : points to multiple rhs entries (ja=0 ; for NRHS=1)
c ia : points back to previously solved dof
c
         ja = (in-1)*MFW
         ia = ja + NFW - 1


c
c recall that Bbuf() contains the elimination accumulated global rhs
c  in this routine b is positioned relative to this dof
c  such that Bbuf(1) = f(m) for rhs1,  Bbuf(2) = f(m) for rhs2,  etc.
c  where m is the dof we are eliminating
c
c for symmetric   : s = f(m)/k(m,m)
c for unsymmetric : s = f(m)
c
         s = Bbuf(in)/f1
         if (iprint.eq.1) then
           write(nfsout,*) 'ELMSOL: in,Bbuf(in) = ',in,Bbuf(in)
           write(nfsout,*) 'ELMSOL: f1 = ',f1
           write(nfsout,*) 'ELMSOL: s = ',s
           call pause
         endif
c
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c for the first pass, there is 1 unknown/1 equation
c  thus neither of these loops gets executed   (ie: x(m) = f(m)/k(m,m))
c
c loop thru the dof beyond this one in the solution front
c --------------------------------------------------------
         do 10 i = idp,NFW



c
c suffle the solution forward 1 slot to make room for this solution
c  if necessary, cuz its possible the new solution just gets stuck at th
c
            Xfrnt(ia+1) = Xfrnt(ia)
c
c compute the contribution of each dof to this solution
c -----------------------------------------------------
c ie: Xfrnt(m) = [ f(m)/k(m,m) - {k(m,m+1)/k(m,m)} * x(m+1)
c                - {k(m,m+2)/k(m,m)} * x(m+2) - ..... ]
c
c matching w/ the equation below:
c  s =  f(m)/k(m,m)
c  Ubuf(iuu) = k(m,m+1)/k(m,m)  (for symmetric)
c  Ubuf(iuu) = k(m,m+1)         (for unsymmetric)
c  Xfrnt(ia) = Xfrnt(m+1)
c
            s = s - Ubuf(iuu)*Xfrnt(ia)
            ia = ia - 1
            iuu = iuu - 1
c
   10    continue
c
         iuu = iuu - 1
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c



         do 30 i = 1,idm
            s = s - Ubuf(iuu)*Xfrnt(ia)
            ia = ia - 1
            iuu = iuu - 1
   30    continue
c
c store the solution
c -------------------
c for unsymmetric, we must also divide thru by k(m,m)
c
         Xfrnt(ja+id) = s/f2
c
   40 continue
c
      return
      end

