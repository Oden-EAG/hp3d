c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: Assembles lhs into the front for unsymmetric matrices
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c
c I   Amdest:  dof destination vectors
c I   Ellhs :  element lhs
c IO  Flhs  :  front assembled lhs
c               note: for unsymmetric, lhs is stored as:
c               a11 a12 a13 ...
c               a21 a22 a23 ...   ==>  [a11,a12,a13,..., a21,a22,a23, ..
c               a31 a32 a33 ...
c                :   :   :
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
      subroutine unsasm (Amdest, Ellhs, Flhs)
c
      use surfsc1
      use surfsc2
c
      implicit none
c
      complex(8) :: Amdest(*), Ellhs(*) ,Flhs(*)
c
      integer, parameter :: ntim=660
      integer :: ntemp(ntim)
c
      integer :: i,j,mdi,mdj,ml,n
c
      n = 1
c
      if (NDOFM .gt. ntim) then
         write(NFSOUT,7771)
7771      format(1x,'THE TEMP ARRAY: ntemp MUST ',
     .           'BE INCREASED IN UNSASM')
         IERR = 1
         return
      endif
c
c loop thru the dof in the element
c
      do 100 i = 1,NDOFM
c     ------------------
c pull the dof destination vector (for this row)
c
         mdi = int(Amdest(i))
c
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
         do 45 j = 1, NDOFM
            ntemp(j) = n + j - 1
  45     continue
c
c ALLIANT directives
cvd$ nodepchk
cvd$ nosync
c ARDENT directives
c$doit IVDEP
cwb <
c
c loop thru the dof in the element
c
        do 50 j = 1,NDOFM
c        -------------
c pull the dof destination vector (for this column)
c
            mdj = int(Amdest(j))
c
cwb >
c determine the position in the front (function loclu)
cwb   (**note: we may wish to inline this code)
c
cwb             ml = loclu(mdi,mdj)
            ml = mdj + (mdi-1)*NFW
cwb <
c
c assemble into the front
cwb >
cwb             Flhs(ml) = Flhs(ml) + Ellhs(n)
cwb             n = n + 1



            Flhs(ml) = Flhs(ml) + Ellhs(ntemp(j))
cwb <
c
   50   continue
cwb >
        n = n + NDOFM
cwb <
c
100   continue
c
c
      end subroutine unsasm
