!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! Calculates the solution for one dof specified by (id)
!
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   Aidx  : position of the dof being solved for within the current fr
!             (more specifically the relative position in u()
!
! I   Ubuf  : the eliminated lhs buffered equations
!              positioned relative to the this equation
!
! I   Bbuf  : the eliminated rhs buffered equations
!              positioned relative to the this equation
!
! IO  Xfrnt : the solution front
cc*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
! LATEST REVISION: Mar 2023
!++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
! NAMING CONVENTIONS:
!     AAAAAAAA    Variables in COMMON & PARAMETERS
!     Aaaaaaaa    Variables as ARGUMENTS
!     aaaaaaaa    LOCAL Variables
!         7xxx    FORMAT Statements
!         9xxx    ERROR Handling
!+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
! Additional notes:
! for a given element the back-substituion will look something like:
!  (note: this implies subsequent call to sub.elmsol
!
! |k11  -  k13 k14  - |       x3 =  f3/k33
! |k21 k22 k23 k24 k25|  ===> x4 = (f4 - k34*x3)/k44
! | -   -  k33  -   - |       x1 = (f1 - k13*x3 - k14*x4)/k11
! | -   -  k34 k44  - |       x5 = (f5 - k53*x3 - k54*x4 - k51*x1)/k55
! |k51  -  k53 k54 k55|       x2 = (f2 - k23*x3 - k24*x4 - k21*x1 - k25*
!
! where u would be stored as:
!  Ubuf = {..., k21 k22 k23 k24 k25, k51 k53 k54 k55, k11 k13 k14, k34 k
!  Bbuf = {..., f2, f5, f1, f4, f3 }
!  Xfrnt= {..., x2, x5, x1, x4, x3 }
!
! note: that it doesnt look like this forever,( ie: this is really a fro
!        some equations remain till the next element, some are 'forgotte
!        they have already been passed to permanent storage, and the 'fo
!        ones become new destinations in Xfrnt() based on Aidx
!+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
!
      subroutine elmsol (Aidx, Ubuf, Bbuf, Xfrnt)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      complex(8) :: Aidx
      complex(8) :: Ubuf(*), Bbuf(*), Xfrnt(*)
!
      complex(8) :: s,f1,f2
      integer    :: i,ia,id,idm,idp,in,iuu,ja
!
      complex(8) :: one = (1.d0,0.d0)
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      id = int(Aidx)
      idm = id - 1
      idp = id + 1
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(nfsout,*) 'ELMSOL: Aidx = ',Aidx
      endif
#endif
!
! because Ubuf contains slightly different info for symmetric versus uns
!  we must take a different path herein
!
!  for symmetric: Ubuf() contains:
!       [k(1,m)/k(m,m), k(2,m)/k(m,m),..., k(m,m), ...,k(nfw,m)/k(m,m)]
!  for unsymmetric: Ubuf() contains:
!       [k(1,m), k(2,m),...k(nfw,m)]
! --------------------------------------------------------------------
! pull the k(m,m) entry from the Ubuf buffer
!-----------------------------------------
      if (ISYM.eq.1 .or. ISYM.eq.4) then
         f1 = Ubuf(id)
         f2 = one
      else
         f1 = one
         f2 = Ubuf(id)
      endif
!
! loop thru and back-substitute for each rhs
! --------------------------------------------
      do 40 in = 1,NRHS
!     -----------------
         iuu = NFW

!
! ja : points to multiple rhs entries (ja=0 ; for NRHS=1)
! ia : points back to previously solved dof
!
         ja = (in-1)*MFW
         ia = ja + NFW - 1

!
! recall that Bbuf() contains the elimination accumulated global rhs
!  in this routine b is positioned relative to this dof
!  such that Bbuf(1) = f(m) for rhs1,  Bbuf(2) = f(m) for rhs2,  etc.
!  where m is the dof we are eliminating
!
! for symmetric   : s = f(m)/k(m,m)
! for unsymmetric : s = f(m)
!
         s = Bbuf(in)/f1
#if HP3D_DEBUG
         if (iprint.eq.1) then
           write(nfsout,*) 'ELMSOL: in,Bbuf(in) = ',in,Bbuf(in)
           write(nfsout,*) 'ELMSOL: f1 = ',f1
           write(nfsout,*) 'ELMSOL: s = ',s
           call pause
         endif
#endif
!
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!
! for the first pass, there is 1 unknown/1 equation
!  thus neither of these loops gets executed   (ie: x(m) = f(m)/k(m,m))
!
! loop thru the dof beyond this one in the solution front
! --------------------------------------------------------
         do 10 i = idp,NFW

!
! suffle the solution forward 1 slot to make room for this solution
!  if necessary, cuz its possible the new solution just gets stuck at th
!
            Xfrnt(ia+1) = Xfrnt(ia)
!
! compute the contribution of each dof to this solution
! -----------------------------------------------------
! ie: Xfrnt(m) = [ f(m)/k(m,m) - {k(m,m+1)/k(m,m)} * x(m+1)
!                - {k(m,m+2)/k(m,m)} * x(m+2) - ..... ]
!
! matching w/ the equation below:
!  s =  f(m)/k(m,m)
!  Ubuf(iuu) = k(m,m+1)/k(m,m)  (for symmetric)
!  Ubuf(iuu) = k(m,m+1)         (for unsymmetric)
!  Xfrnt(ia) = Xfrnt(m+1)
!
            s = s - Ubuf(iuu)*Xfrnt(ia)
            ia = ia - 1
            iuu = iuu - 1
!
   10    continue
!
         iuu = iuu - 1
!
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!

         do 30 i = 1,idm
            s = s - Ubuf(iuu)*Xfrnt(ia)
            ia = ia - 1
            iuu = iuu - 1
   30    continue
!
! store the solution
! -------------------
! for unsymmetric, we must also divide thru by k(m,m)
!
         Xfrnt(ja+id) = s/f2
!
   40 continue
!
!
      end subroutine elmsol
