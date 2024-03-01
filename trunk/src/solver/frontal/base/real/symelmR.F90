      
!
!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: Elimination of one equation (id) for symmetric matricies
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
!
! I   Iel   : element number
!
! I   Aidx  : transfer loc for active dof in the front after eliminating
!             this is the equation to be eliminated
!
! IO  Flhs  : the assembled front lhs
!             note: for symmetric, lhs is stored as:  *** by column ***
!               a11 a12 a13 ...
!                   a22 a23 ...      ==>   [a11, a12,a22, a13,a23,a33, .
!                       a33 ...
!                            :
!
! IO  Ubuf     : the buffer to hold eliminated lhs equations
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
   subroutine symelm (Iel, Aidx, Flhs, Ubuf)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      integer :: Iel
      real(8) :: Aidx
      real(8) :: Flhs(*), Ubuf(*)
!
      real(8) :: pivot,s
      integer :: i,id,idm,idp,j,k,m,mp,n,nn
!
      real(8), parameter :: sml   = 1.d-30
      real(8), parameter :: dzero = 0.0d0
      real(8), parameter :: sml2  = 1.d-15
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
!
      if (iprint.eq.1) then
        write(*,*) 'SYMELM: Iel,Aidx,Flhs,Ubuf = ',Iel,Aidx
        write(*,7001) (Flhs(iii),iii=1,10)
        write(*,7001) (Ubuf(iii),iii=1,10)
 7001   format(1x,5e12.5)
        call pause
      endif
#endif
!
! set up pointers into Flhs
!
      id = nint(Aidx)
      mp = (id*(id+1))/2
      idm = id - 1
      idp = id + 1
      m = mp - id + 1
      k = 1
!
! pull the diagonal entry
!
      pivot = Flhs(mp)
!
#if HP3D_DEBUG
      if(IPRPIV .eq. 1 .and. iprint .eq. 1) then
         write(NFSOUT,1000) Iel,nfw,id,pivot
 1000    format(5x,'IEL,NFW,ID,PIVOT',3i5,1pe11.3)
      endif
#endif
!
      Ubuf(id) = pivot
!
! check for zero pivots
!
      if (dabs(pivot) .le. sml) then
         IERR=2
         return
      endif
!
      if(pivot .lt. dzero) NNEGP = NNEGP+1
      if(pivot .gt. dzero) NPOSP = NPOSP+1
!
! perform the lhs elimination:
!-----------------------------
!   k(i,j)' = k(i,j) - [k(i,m)/k(m,m)]*k(m,j)
!
!  where m is the row being eliminated
!
! loop thru the part of Flhs above this dof and perform:
! ------------------------------------------------------
! into u() we place:  [k(1,m)/k(m,m), k(2,m)/k(m,m),...k(m-1,m)/k(m,m),
!
      do 20 i = 1,idm
!     ---------------
! pick up k(i,m)
!
         s = Flhs(m)
!
! compute k(i,m)/k(m,m)
!
         Ubuf(i) = s/pivot
!wb >
! dont operate on zero elements
! -----------------------------
         if (dabs(s) .le. sml2) then
            m = m + 1
            k = k + i
            go to 20
          endif
!wb <
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
! eliminate
! ----------
!   k(i,j)' = k(i,j) - [k(i,m)/k(m,m)]*k(m,j)      note: k(i,m) = k(m,i)
!
         do 10 j=1,i
            Flhs(k) = Flhs(k) - s*Ubuf(j)
            k = k + 1
   10    continue
!
         m = m + 1
   20 continue
!
!-----------------------------------------------------------------------
!
! loop thru the dof below the elimination dof in the front
!----------------------------------------------------------
! now into u() we place: [... k(m,m), k(m,m+1)/k(m,m), k(m,m+2)/k(m,m),.
!
      m = mp
      k = 0
      do 80 i = idp,NFW
!     -----------------
         nn = m - id
         m = m + id + k
         n = m - id
!
! pick up k(m,i)
!
         s = Flhs(m)
!
! compute k(m,i)/k(m,m)
!
         Ubuf(i) = s/pivot
!wb >
! dont operate on zero elements, just shuffle entries forward
! ------------------------------------------------------------
         if (dabs(s) .le. sml2) then
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
            do 30 j = 1,idm
               Flhs(nn+j) = Flhs(n+j)
   30       continue
            nn = nn - 1
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
            do 35 j = idp,i
               Flhs(nn+j) = Flhs(n+j)
   35       continue
!
            k = k + 1
            go to 80
          endif
!wb <
!
! eliminate
! ---------
!   k(i,j)' = k(i,j) - [k(i,m)/k(m,m)]*k(m,j)
!
! **note: we are suffling the equations forward as we perform the calcs
!         the suffling occurs as:  say we eliminate column 3
!
!        a11 a12 a13 a14 a15       a11 a12  *  a14 a15      a11 a12 a14
!            a22 a23 a24 a25           a22  *  a24 a25          a22 a24
!                a33 a34 a35  ==>           x   *   *  ==>          a44
!                    a44 a45                   a44 a45
!                        a55                       a55
!
!  note: the * entries are the k(m,i) signified above
!
!  ie: k14' = k14 - [k13/k33]*k34
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
         do 40 j = 1,idm
            Flhs(nn+j) = Flhs(n+j) - s*Ubuf(j)
   40    continue
!
         nn = nn - 1
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
         do 60 j = idp,i
            Flhs(nn+j) = Flhs(n+j) - s*Ubuf(j)
   60    continue
!
         k = k + 1
!
   80 continue
!  ------------
!
   end subroutine symelm
