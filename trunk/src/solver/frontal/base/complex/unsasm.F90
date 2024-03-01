!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: Assembles lhs into the front for unsymmetric matrices
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
!
! I   Amdest:  dof destination vectors
! I   Ellhs :  element lhs
! IO  Flhs  :  front assembled lhs
!               note: for unsymmetric, lhs is stored as:
!               a11 a12 a13 ...
!               a21 a22 a23 ...   ==>  [a11,a12,a13,..., a21,a22,a23, ..
!               a31 a32 a33 ...
!                :   :   :
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
      subroutine unsasm (Amdest, Ellhs, Flhs)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      complex(8) :: Amdest(*), Ellhs(*) ,Flhs(*)
!
      integer, parameter :: ntim=660
      integer :: ntemp(ntim)
!
      integer :: i,j,mdi,mdj,ml,n
!
      n = 1
!
      if (NDOFM .gt. ntim) then
         write(NFSOUT,7771)
7771      format(1x,'THE TEMP ARRAY: ntemp MUST ', &
                 'BE INCREASED IN UNSASM')
         IERR = 1
         return
      endif
!
! loop thru the dof in the element
!
      do 100 i = 1,NDOFM
!     ------------------
! pull the dof destination vector (for this row)
!
         mdi = int(Amdest(i))
!
! ALLIANT directives
!cvd$ select (vector)
!cvd$ nodepchk
!cvd$ nosync
! ARDENT directives
!$doit VBEST
!$doit IVDEP
!
!
! to remove vector dependencies, we fill a temporary vector with values
!
         do 45 j = 1, NDOFM
            ntemp(j) = n + j - 1
  45     continue
!
! ALLIANT directives
!cvd$ nodepchk
!cvd$ nosync
! ARDENT directives
!$doit IVDEP
!cwb <
!
! loop thru the dof in the element
!
        do 50 j = 1,NDOFM
!        -------------
! pull the dof destination vector (for this column)
!
            mdj = int(Amdest(j))
!
!cwb >
! determine the position in the front (function loclu)
!cwb   (**note: we may wish to inline this code)
!
!cwb             ml = loclu(mdi,mdj)
            ml = mdj + (mdi-1)*NFW
!cwb <
!
! assemble into the front
!cwb >
!cwb             Flhs(ml) = Flhs(ml) + Ellhs(n)
!cwb             n = n + 1



            Flhs(ml) = Flhs(ml) + Ellhs(ntemp(j))
!cwb <
!
   50   continue
!cwb >
        n = n + NDOFM
!cwb <
!
100   continue
!
!
      end subroutine unsasm