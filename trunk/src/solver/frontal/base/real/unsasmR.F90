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
      real(8) :: Amdest(*), Ellhs(*), Flhs(*)
      integer :: ntemp(324)
!
      integer :: n,i,j,mdi,mdj,ml
!
      n = 1
!
      if (NDOFM .gt. 324) then
         write(NFSOUT,7771)
7771      format(1x,'THE TEMP ARRAY: ntemp MUST ', &
                 'BE INCREASED IN UNSASM')
         IERR = 1
         return
      endif
!
! loop thru the dof in the element
!
      do i = 1,NDOFM
!     ------------------
! pull the dof destination vector (for this row)
!
         mdi = int(Amdest(i))
!wb >
! ALLIANT directives
!vd$ select (vector)
!vd$ nodepchk
!vd$ nosync
! ARDENT directives
!$doit VBEST
!$doit IVDEP
!
!
! to remove vector dependencies, we fill a temporary vector with values
!
         do j = 1, NDOFM
            ntemp(j) = n + j - 1
         enddo
!
! ALLIANT directives
!vd$ nodepchk
!vd$ nosync
! ARDENT directives
!$doit IVDEP
!wb <
!
! loop thru the dof in the element
!
        do j = 1,NDOFM
!        -------------
! pull the dof destination vector (for this column)
!
            mdj = int(Amdest(j))
!
!wb >
! determine the position in the front (function loclu)
!wb   (**note: we may wish to inline this code)
!
!wb             ml = loclu(mdi,mdj)
            ml = mdj + (mdi-1)*NFW
!wb <
!
! assemble into the front
!wb >
!wb             Flhs(ml) = Flhs(ml) + Ellhs(n)
!wb             n = n + 1
            Flhs(ml) = Flhs(ml) + Ellhs(ntemp(j))
!wb <
!
        enddo
!wb >
        n = n + NDOFM
!wb <
!
      enddo
!
!
   end subroutine unsasm
