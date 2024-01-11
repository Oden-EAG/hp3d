      
c
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: Elimination of one equation (id) for symmetric matricies
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c
c I   Iel   : element number
c
c I   Aidx  : transfer loc for active dof in the front after eliminating
c             this is the equation to be eliminated
c
c IO  Flhs  : the assembled front lhs
c             note: for symmetric, lhs is stored as:  *** by column ***
c               a11 a12 a13 ...
c                   a22 a23 ...      ==>   [a11, a12,a22, a13,a23,a33, .
c                       a33 ...
c                            :
c
c IO  Ubuf     : the buffer to hold eliminated lhs equations
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
      subroutine symelm (Iel, Aidx, Flhs, Ubuf)
c
      use surfsc1
      use surfsc2
c
      implicit none
c
      integer :: Iel
      real(8) :: Aidx
      real(8) :: Flhs(*), Ubuf(*)
c
      real(8) :: pivot,s
      integer :: i,id,idm,idp,j,k,m,mp,n,nn
c
      real(8), parameter :: sml   = 1.d-30
      real(8), parameter :: dzero = 0.0d0
      real(8), parameter :: sml2  = 1.d-15
c
#if DEBUG_MODE
      integer :: iprint
      iprint=0
c
      if (iprint.eq.1) then
        write(*,*) 'SYMELM: Iel,Aidx,Flhs,Ubuf = ',Iel,Aidx
        write(*,7001) (Flhs(iii),iii=1,10)
        write(*,7001) (Ubuf(iii),iii=1,10)
 7001   format(1x,5e12.5)
        call pause
      endif
#endif
c
c set up pointers into Flhs
c
      id = nint(Aidx)
      mp = (id*(id+1))/2
      idm = id - 1
      idp = id + 1
      m = mp - id + 1
      k = 1
c
c pull the diagonal entry
c
      pivot = Flhs(mp)
c
#if DEBUG_MODE
      if(IPRPIV .eq. 1 .and. iprint .eq. 1) then
         write(NFSOUT,1000) Iel,nfw,id,pivot
 1000    format(5x,'IEL,NFW,ID,PIVOT',3i5,1pe11.3)
      endif
#endif
c
      Ubuf(id) = pivot
c
c check for zero pivots
c
      if (dabs(pivot) .le. sml) then
         IERR=2
         return
      endif
c
      if(pivot .lt. dzero) NNEGP = NNEGP+1
      if(pivot .gt. dzero) NPOSP = NPOSP+1
c
c perform the lhs elimination:
c-----------------------------
c   k(i,j)' = k(i,j) - [k(i,m)/k(m,m)]*k(m,j)
c
c  where m is the row being eliminated
c
c loop thru the part of Flhs above this dof and perform:
c ------------------------------------------------------
c into u() we place:  [k(1,m)/k(m,m), k(2,m)/k(m,m),...k(m-1,m)/k(m,m),
c
      do 20 i = 1,idm
c     ---------------
c pick up k(i,m)
c
         s = Flhs(m)
c
c compute k(i,m)/k(m,m)
c
         Ubuf(i) = s/pivot
cwb >
c dont operate on zero elements
c -----------------------------
         if (dabs(s) .le. sml2) then
            m = m + 1
            k = k + i
            go to 20
          endif
cwb <
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c eliminate
c ----------
c   k(i,j)' = k(i,j) - [k(i,m)/k(m,m)]*k(m,j)      note: k(i,m) = k(m,i)
c
         do 10 j=1,i
            Flhs(k) = Flhs(k) - s*Ubuf(j)
            k = k + 1
   10    continue
c
         m = m + 1
   20 continue
c
c-----------------------------------------------------------------------
c
c loop thru the dof below the elimination dof in the front
c----------------------------------------------------------
c now into u() we place: [... k(m,m), k(m,m+1)/k(m,m), k(m,m+2)/k(m,m),.
c
      m = mp
      k = 0
      do 80 i = idp,NFW
c     -----------------
         nn = m - id
         m = m + id + k
         n = m - id
c
c pick up k(m,i)
c
         s = Flhs(m)
c
c compute k(m,i)/k(m,m)
c
         Ubuf(i) = s/pivot
cwb >
c dont operate on zero elements, just shuffle entries forward
c ------------------------------------------------------------
         if (dabs(s) .le. sml2) then
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
            do 30 j = 1,idm
               Flhs(nn+j) = Flhs(n+j)
   30       continue
            nn = nn - 1
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
            do 35 j = idp,i
               Flhs(nn+j) = Flhs(n+j)
   35       continue
c
            k = k + 1
            go to 80
          endif
cwb <
c
c eliminate
c ---------
c   k(i,j)' = k(i,j) - [k(i,m)/k(m,m)]*k(m,j)
c
c **note: we are suffling the equations forward as we perform the calcs
c         the suffling occurs as:  say we eliminate column 3
c
c        a11 a12 a13 a14 a15       a11 a12  *  a14 a15      a11 a12 a14
c            a22 a23 a24 a25           a22  *  a24 a25          a22 a24
c                a33 a34 a35  ==>           x   *   *  ==>          a44
c                    a44 a45                   a44 a45
c                        a55                       a55
c
c  note: the * entries are the k(m,i) signified above
c
c  ie: k14' = k14 - [k13/k33]*k34
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
         do 40 j = 1,idm
            Flhs(nn+j) = Flhs(n+j) - s*Ubuf(j)
   40    continue
c
         nn = nn - 1
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
         do 60 j = idp,i
            Flhs(nn+j) = Flhs(n+j) - s*Ubuf(j)
   60    continue
c
         k = k + 1
c
   80 continue
c  ------------
c
      end subroutine symelm
