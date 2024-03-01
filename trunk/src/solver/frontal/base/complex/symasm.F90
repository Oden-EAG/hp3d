!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: Assembles the lhs into the front for symmetric matrices
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   Amdest:  dof destination vectors
! I   Ellhs :  element lhs
! IO  Flhs  :  front assembled lhs
!               note: for symmetric, lhs is stored as:  *** by column **
!               a11 a12 a13 ...
!                   a22 a23 ...      ==>   [a11, a12,a22, a13,a23,a33, .
!                       a33 ...
!                            :
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
      subroutine symasm (Amdest, Ellhs, Flhs)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      complex(8) :: Amdest(*), Ellhs(*), Flhs(*)
!
      integer, parameter :: ntim=324
      integer :: ntemp(ntim)
!
      integer :: i,j,k,l,mi,mj,mk,n
!
      n = 1
!
      if (NDOFM .gt. ntim) then
         write(NFSOUT,7771)
7771      format(1x,'THE TEMP ARRAY: ntemp MUST ', &
                 'BE INCREASED IN SYMASM')
         IERR = 1
         return
      endif
!
!
! loop thru the dof in the element
!
      do 100 i = 1,NDOFM
!     ------------------
! pull the dof destination vector (for this row)
!
          mi = int(Amdest(i))
!cwb >
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
         do 45 j = 1, i
            ntemp(j) = n + j - 1
  45     continue
!
! ALLIANT directives
!cvd$ nodepchk
!cvd$ nosync
! ARDENT directives
!$doit IVDEP
!cwb <
! becuz its symmetric loop thru lower triangular dof for this element
!
         do 50 j = 1, i
!        -------------
! pull the dof destination vector (for this column)
!
            mj = int(Amdest(j))
!cwb >
! determine the position in the front (function locls)
!cwb   (**note: we may wish to inline this code)
!
!cwb             mk = locls(mi,mj)
            k = max0(mi,mj)
            l = min0(mi,mj)
!
            mk = k*(k-1)/2 + l
!cwb <
!
! assemble into the front
!cwb >
!cwb             Flhs(mk) = Flhs(mk) + Ellhs(n)
!cwb             n = n + 1
            Flhs(mk) = Flhs(mk) + Ellhs(ntemp(j))
!cwb <
!
 50      continue
!cwb >
         n = n + i
!cwb <
!
100   continue
!
!
      end subroutine symasm
