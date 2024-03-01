!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: Elimination of one equation (id) for unsymmetric matrices
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   Iel   : element number
!
! I   Aidx  : transfer loc for active dof in the front after eliminating
!             this is the equation to be eliminated
!
! IO  Flhs  : the assembled front lhs
!             note: for unsymmetric, lhs is stored as:
!               a11 a12 a13 ...
!               a21 a22 a23 ...   ==>   [a11,a12,a13,..., a21,a22,a23, .
!               a31 a32 a33 ...
!                :   :   :
!
! IO  Ubuf   : the buffer to hold eliminated lhs equations
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
      subroutine unselm (Iel, Aidx, Flhs, Ubuf)
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
      integer :: i,id,idm,idp,j,k,kk,m,mp
!
      real(8), parameter :: sml   = 1.d-30
      real(8), parameter :: dzero = 0.0d0
      real(8), parameter :: sml2  = 1.d-15
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
! set up pointers into Flhs
!
      id = int(Aidx)
      idm = id - 1
      idp = id + 1
!
! set the offset for row manipulation
!   (ie: each equation has KFW entries, so id begins at (id-1)*KFW)
!
      k = idm*KFW
      mp = k + id
!
! pull the diagonal entry
!
      pivot = Flhs(mp)
!
#if HP3D_DEBUG
      if(IPRPIV .eq. 1 .and. iprint .eq. 1) then
         write(NFSOUT,7000)  Iel,NFW,id,pivot
 7000    format(5x,'IEL,NFW,ID,PIVOT',3i5,1pe11.3)
      endif
#endif
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
!cvd$ select (vector)
!
!*** for the unsymmetric case we copy the assembled Flhs row into Ubuf()
!====================================================================
! ie: into u() we place:  [k(1,m), k(2,m),...k(NFW,m)]
!     (for the symmetric case u is divided by the pivot)
!    *note: these are of course the k(i,j)' from previous eliminations
!
      do 5 i = 1,NFW
         Ubuf(i) = Flhs(k+i)
    5 continue
!
!-----------------------------------------------------------------------
! perform the lhs elimination:
!-----------------------------
!   k(i,j)' = k(i,j) - [k(i,m)/k(m,m)]*k(m,j)
!
!  where m is the row being eliminated
!
! **note: we are suffling the equations forward as we perform the calcs
!         the suffling occurs as:  say we eliminate column 3
!
!  a11 a12 a13 a14 a15     1  :   a11 a12  #  a14 a15                a11
!  a21 a22 a23 a24 a25     idm:   a21 a22  #  a24 a25     b1 # b2    a21
!  a31 a32 a33 a34 a35  ==>id :    *   *   x   *   * ==>   * x * ==> a31
!  a41 a42 a43 a44 a45     idp:   a41 a42  #  a44 a45     b3 # b4    a51
!  a51 a52 a53 a54 a55     NFW:   a51 a52  #  a54 a55     =======
!
!  note: the * entries are the k(m,i) signified above
!
!  ie: k14' = k14 - [k13/k33]*k34
!
! loop thru the part of Flhs above this dof and perform:
!
      k = 0
      kk = 0
!
      do 30 i = 1,idm
!     ----------------
! *note: (id+k) is an appropriate offset to pick up the column entries (
!        (j+k)  is an appropriate offset to pick up the row entries (a11
!        and Ubuf(j) represents * above
!
         s = Flhs(id+k)/pivot
!cwb >
! dont operate on zero elements
! -----------------------------
         if (dabs(s) .le. sml2) then
            m = k + idm
         else
!cwb <
! elimination for block b1 (as shown in the figure above)
! ------------------------
!cvd$ select (vector)
!
            do 10 j = 1,idm
               m = j + k
               Flhs(m) = Flhs(m) - s*Ubuf(j)
   10       continue
         endif
!
         m = k - 1
!
!cwb >
! dont operate on zero elements
! -----------------------------
         if (dabs(s) .le. sml2) then
!
!cvd$ select (vector)
!
            do 15 j = idp,NFW
               Flhs(j+m) = Flhs(j+k)
   15       continue
         else
!cwb <
! elimination for block b2 (as shown in the figure above)
! ------------------------
!cvd$ select (vector)
!
            do 20 j = idp,NFW
               Flhs(j+m) = Flhs(j+k) - s*Ubuf(j)
   20       continue
         endif
!
! store also the k(i,m)/k(m,m) terms, for possible dumping to diskfile
!  (to the l buffer) for unsymmetric w/ resolution
!  and for use in elmrhs, to eliminate the rhs
!
         k = k + KFW
         Flhs(k-KFW+NFW) = s
!
   30 continue
!
!-----------------------------------------------------------------------
! loop thru the dof below the elimination dof in the front
!
      k= k + KFW
      do 70 i = idp,NFW
!    -------------------
         s = Flhs(id+k)/pivot
         m = k - KFW
!cwb >
! dont operate on zero elements
! -----------------------------
         if (dabs(s) .le. sml2) then
!
!cvd$ select (vector)
!
            do 40 j = 1,idm
              Flhs(j+m) = Flhs(k+j)
   40       continue
         else
!cwb <
!
! elimination for block b3 (as shown in the figure above)
! ------------------------
!cvd$ select (vector)
!
            do 50 j = 1,idm
              Flhs(j+m) = Flhs(k+j) - s*Ubuf(j)
   50       continue
         endif
!
         m = m - 1
!
!cwb >
! dont operate on zero elements
! -----------------------------
         if (dabs(s) .le. sml2) then
!
!cvd$ select (vector)
!
            do 55 j = idp,NFW
               Flhs(m+j) = Flhs(k+j)
   55       continue
         else
!cwb <
!
! elimination for block b4 (as shown in the figure above)
! ------------------------
!cvd$ select (vector)
!
            do 60 j = idp,NFW
               Flhs(m+j) = Flhs(k+j) - s*Ubuf(j)
   60       continue
         endif
!
! store also the k(i,m)/k(m,m) terms, for possible dumping to diskfile
!  (to the l buffer) for unsymmetric w/ resolution
!  and for use in elmrhs, to eliminate the rhs
!
         Flhs(k-KFW+NFW) = s
         k = k + KFW
!
   70 continue
!
!
      end subroutine unselm
