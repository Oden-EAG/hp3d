!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! Elimination of rhs's for equation (id)
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
!
! I   Aidx   transfer loc for active dof in the front after eliminating
!             this is the equation to be eliminated
!
! I   Inc    for unsymmetric this is the Increment to shift for picking
!             the correct entries in Ubuf()
!
! IO  Frhs   the assembled front rhs
!
! IO  Ubuf   the buffer to hold eliminated lhs equations
! IO  Bbuf   the buffer to hold eliminated rhs equations
!
!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
! LATEST REVISION: Mar 2023
!++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
! NAMING CONVENTIONS:
!     AAAAAAAA    Variables in COMMON & PARAMETERS
!     Aaaaaaaa    Variables as ARGUMENTS
!     aaaaaaaa    LOCAL Variables
!         7xxx    FORMAT Statements
!         9xxx    ERROR Handling
!+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
!
      subroutine elmrhs (Aid, Inc, Frhs, Ubuf, Bbuf)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      real(8) :: Aid
      integer :: Inc
      real(8) :: Frhs(*), Ubuf(*), Bbuf(*)
!
      real(8) :: s
      integer :: in,id,idm,idp,ii,ii1,ii2,im,iuu
!
      real(8), parameter :: sml2 = 1.d-15
!
      id = int(Aid)
      idm = id - 1
      idp = id + 1
      im = 0
!
! perform the rhs elimination:
!-----------------------------
!   f(i)' = f(i) - [k(i,m)/k(m,m)]*f(m)
!
!  where m is the row being eliminated
!
!
! loop over the number of rhs in the model
!
      do 60 in = 1,NRHS
!     -----------------
!
         iuu = 1
!
! pull f(m)
!
         s = Frhs(im+id)
!
! ***note: we load the rhs elimination buffer with the last  f(i)'
!           (ie: f(i)" = (f(i)'-[k(i,m)/k(m,m)]*f(m))....)
!
         Bbuf(in) = s
!        ----------
!cwb >
! dont operate on zero elements
! -----------------------------
         if (dabs(s) .le. sml2) then
            iuu = iuu + idm*Inc
            go to 25
         endif
!cwb <
!
         ii1 = im + 1
         ii2 = im + idm
!
! loop over the front dof, up to the dof we are eliminating
! note: u() contains already the entries [k(i,m)/k(m,m)]
! and eliminate
!     ---------
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!
!
         do 20 ii = ii1,ii2
            Frhs(ii) = Frhs(ii) - s*Ubuf(iuu)
            iuu = iuu + Inc
   20    continue
!
   25    if (ISYM.eq.1 .or. ISYM.eq.4) iuu = iuu + 1
!
!cwb          if (idp .gt. NFW) go to 50
!
! loop thru the dof below the elimination dof in the front
!----------------------------------------------------------
! **note: we are suffling the equations forward as we perform the calcs
!
        if (idp .le. NFW) then
!
            ii1 = im + idp
            ii2 = im + NFW
!
!cwb >
! dont operate on zero elements
! -----------------------------
            if (dabs(s) .le. sml2) then
!
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!
!
               do 27 ii = ii1,ii2
                  Frhs(ii-1) = Frhs(ii)
   27          continue
!
               iuu = iuu + (NFW - idp + 1)*Inc
            else
!cwb <
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!
!
               do 30 ii = ii1,ii2
                  Frhs(ii-1) = Frhs(ii) - s*Ubuf(iuu)
                  iuu = iuu + Inc
   30          continue
            endif
         endif
!
         im = im + MFW
   60 continue
!
!
      end subroutine elmrhs
