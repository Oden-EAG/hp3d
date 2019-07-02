      subroutine symasm (Amdest, Ellhs, Flhs)
c
#include "syscom.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: Assembles the lhs into the front for symmetric matrices
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Amdest:  dof destination vectors
c I   Ellhs :  element lhs
c IO  Flhs  :  front assembled lhs
c               note: for symmetric, lhs is stored as:  *** by column **
c               a11 a12 a13 ...
c                   a22 a23 ...      ==>   [a11, a12,a22, a13,a23,a33, .
c                       a33 ...
c                            :
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
ctim      dimension Amdest(*), Ellhs(*), Flhs(*), ntemp(100)
ccccc dimension Amdest(*), Ellhs(*), Flhs(*), ntemp(324)
      dimension Amdest(*), Ellhs(*), Flhs(*), ntemp(10000)
c
      n = 1
c
ctim      if (NDOFM .gt. 100) then
cccc      if (NDOFM .gt. 324) then
      if (NDOFM .gt. 10000) then
         write(NFSOUT,7771)
7771      format(1x,'THE TEMP ARRAY: ntemp MUST ',
     .           'BE INCREASED IN SYMASM')
         IERR = 1
         return
      endif



c      write(*,*) 'in symasm ellhs = '
c      kk=0
c      do 44 ii=1,NDOFM
c        do 43 jj=1,ii
c          kk=kk+1
c          write(*,*) Ellhs(kk)
c   43   continue
c        pause
c   44 continue


      
c
c loop thru the dof in the element
c
      do 100 i = 1,NDOFM
c     ------------------
c pull the dof destination vector (for this row)
c
          mi = Amdest(i)
cwb >
c ALLIANT directives
cvd$ select (vector)
cvd$ nodepchk
cvd$ nosync
c ARDENT directives
c$doit VBEST
c$doit IVDEP
c
c
c to remove vector dependencies, we fill a temporary vector with values
c
         do 45 j = 1, i
            ntemp(j) = n + j - 1
  45     continue
c
c ALLIANT directives
cvd$ nodepchk
cvd$ nosync
c ARDENT directives
c$doit IVDEP
cwb <
c becuz its symmetric loop thru lower triangular dof for this element
c
         do 50 j = 1, i
c        -------------
c pull the dof destination vector (for this column)
c
            mj = Amdest(j)
cwb >
c determine the position in the front (function locls)
cwb   (**note: we may wish to inline this code)
c
cwb             mk = locls(mi,mj)
            k = max0(mi,mj)
            l = min0(mi,mj)
c
            mk = k*(k-1)/2 + l
cwb <
c
c assemble into the front
cwb >
cwb             Flhs(mk) = Flhs(mk) + Ellhs(n)
cwb             n = n + 1
            Flhs(mk) = Flhs(mk) + Ellhs(ntemp(j))
cwb <
c
 50      continue
cwb >
         n = n + i
cwb <
c
100   continue
c
      return
      end
