      subroutine unselm (Iel, Aidx, Flhs ,Ubuf)
c
#include "syscomz.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: Elimination of one equation (id) for unsymmetric matrices
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Iel   : element number
c
c I   Aidx  : transfer loc for active dof in the front after eliminating
c             this is the equation to be eliminated
c
c IO  Flhs  : the assembled front lhs
c             note: for unsymmetric, lhs is stored as:
c               a11 a12 a13 ...
c               a21 a22 a23 ...   ==>   [a11,a12,a13,..., a21,a22,a23, .
c               a31 a32 a33 ...
c                :   :   :
c
c IO  Ubuf   : the buffer to hold eliminated lhs equations
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
c
      dimension Flhs(*), Ubuf(*)
c
      data sml,zero /1.e-30, 0.0/
cwb >
      data sml2 /1.e-15/
cwb <

      iprint=0


c
c set up pointers into Flhs
c
      id = Aidx
      idm = id - 1
      idp = id + 1
c
c set the offset for row manipulation
c   (ie: each equation has KFW entries, so id begins at (id-1)*KFW)
c
      k = idm*KFW
      mp = k + id
      if (iprint.eq.1) then
        write(*,7001) Aidx,id,idm,idp,KFW,k,mp
        write(*,7002) (Flhs(mp+ii),ii=-3,3)
 7001   format('unselm: Aidx,id,idm,idp,KFW,k,mp = ', 2e12.5,2x,6i4)
 7002   format('unselm: Flhs(mp+ii),ii=-3,3 = ',7(2e12.5,2x))
        call pause
      endif
c
c pull the diagonal entry
c
      pivot = Flhs(mp)
c
      iprint = 0
      if(IPRPIV .eq. 1 .and. iprint .eq. 1) then
         write(NFSOUT,7000)  Iel,NFW,id,pivot
 7000    format(5x,'IEL,NFW,ID,PIVOT',3i5,2e25.15)
      endif
c
c check for zero pivots
c
cwr10.07.99
ccccc if (Cdabs(pivot) .le. CDABS(sml)) then
      if (  abs(pivot) .le.   ABS(sml)) then
         IERR=2
         return
      endif
c
cwr10.07.99
ccccc if(CDABS(pivot) .lt. CDABS(zero)) NNEGP = NNEGP+1
ccccc if(CDABS(pivot) .gt. CDABS(zero)) NPOSP = NPOSP+1
      if(  ABS(pivot) .lt.   ABS(zero)) NNEGP = NNEGP+1
      if(  ABS(pivot) .gt.   Abs(zero)) NPOSP = NPOSP+1
c
cvd$ select (vector)
c
c*** for the unsymmetric case we copy the assembled Flhs row into Ubuf()
c====================================================================
c ie: into u() we place:  [k(1,m), k(2,m),...k(NFW,m)]
c     (for the symmetric case u is divided by the pivot)
c    *note: these are of course the k(i,j)' from previous eliminations
c
      do 5 i = 1,NFW
         Ubuf(i) = Flhs(k+i)
    5 continue
c
c-----------------------------------------------------------------------
c perform the lhs elimination:
c-----------------------------
c   k(i,j)' = k(i,j) - [k(i,m)/k(m,m)]*k(m,j)
c
c  where m is the row being eliminated
c
c **note: we are suffling the equations forward as we perform the calcs
c         the suffling occurs as:  say we eliminate column 3
c
c  a11 a12 a13 a14 a15     1  :   a11 a12  #  a14 a15                a11
c  a21 a22 a23 a24 a25     idm:   a21 a22  #  a24 a25     b1 # b2    a21
c  a31 a32 a33 a34 a35  ==>id :    *   *   x   *   * ==>   * x * ==> a31
c  a41 a42 a43 a44 a45     idp:   a41 a42  #  a44 a45     b3 # b4    a51
c  a51 a52 a53 a54 a55     NFW:   a51 a52  #  a54 a55     =======
c
c  note: the * entries are the k(m,i) signified above
c
c  ie: k14' = k14 - [k13/k33]*k34
c
c loop thru the part of Flhs above this dof and perform:
c
      k = 0
      kk = 0
c
      do 30 i = 1,idm
c     ----------------
c *note: (id+k) is an appropriate offset to pick up the column entries (
c        (j+k)  is an appropriate offset to pick up the row entries (a11
c        and Ubuf(j) represents * above
c
         s = Flhs(id+k)/pivot
cwb >
c dont operate on zero elements
c -----------------------------
cwr10.07.99
ccccc    if (Cdabs(s) .le. CDABS(sml2)) then
         if (  abs(s) .le.   ABS(sml2)) then
            m = k + idm
         else
cwb <
c elimination for block b1 (as shown in the figure above)
c ------------------------
cvd$ select (vector)
c
            do 10 j = 1,idm
               m = j + k
               Flhs(m) = Flhs(m) - s*Ubuf(j)
   10       continue
         endif
c
         m = k - 1
c
cwb >
c dont operate on zero elements
c -----------------------------
cwr10.07.99
ccccc    if (Cdabs(s) .le. CDABS(sml2)) then
         if (  abs(s) .le.   ABS(sml2)) then
c
cvd$ select (vector)
c
            do 15 j = idp,NFW
               Flhs(j+m) = Flhs(j+k)
   15       continue
         else
cwb <
c elimination for block b2 (as shown in the figure above)
c ------------------------
cvd$ select (vector)
c
            do 20 j = idp,NFW
               Flhs(j+m) = Flhs(j+k) - s*Ubuf(j)
   20       continue
         endif
c
c store also the k(i,m)/k(m,m) terms, for possible dumping to diskfile
c  (to the l buffer) for unsymmetric w/ resolution
c  and for use in elmrhs, to eliminate the rhs
c
         k = k + KFW
         Flhs(k-KFW+NFW) = s
c
   30 continue
c
c-----------------------------------------------------------------------
c loop thru the dof below the elimination dof in the front
c
      k= k + KFW
      do 70 i = idp,NFW
c    -------------------
         s = Flhs(id+k)/pivot
         m = k - KFW
cwb >
c dont operate on zero elements
c -----------------------------
cwr10.07.99
ccccc    if (Cdabs(s) .le. CDABS(sml2)) then
         if (  abs(s) .le.   ABS(sml2)) then
c
cvd$ select (vector)
c
            do 40 j = 1,idm
              Flhs(j+m) = Flhs(k+j)
   40       continue
         else
cwb <
c
c elimination for block b3 (as shown in the figure above)
c ------------------------
cvd$ select (vector)
c
            do 50 j = 1,idm
              Flhs(j+m) = Flhs(k+j) - s*Ubuf(j)
   50       continue
         endif
c
         m = m - 1
c
cwb >
c dont operate on zero elements
c -----------------------------
cwr10.07.99
ccccc    if (Cdabs(s) .le. CDABS(sml2)) then
         if (  abs(s) .le.   ABS(sml2)) then
c
cvd$ select (vector)
c
            do 55 j = idp,NFW
               Flhs(m+j) = Flhs(k+j)
   55       continue
         else
cwb <
c
c elimination for block b4 (as shown in the figure above)
c ------------------------
cvd$ select (vector)
c
            do 60 j = idp,NFW
               Flhs(m+j) = Flhs(k+j) - s*Ubuf(j)
   60       continue
         endif
c
c store also the k(i,m)/k(m,m) terms, for possible dumping to diskfile
c  (to the l buffer) for unsymmetric w/ resolution
c  and for use in elmrhs, to eliminate the rhs
c
         Flhs(k-KFW+NFW) = s
         k = k + KFW
c
   70 continue
c
      return
      end
