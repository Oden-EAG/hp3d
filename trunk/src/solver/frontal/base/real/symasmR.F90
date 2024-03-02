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
      real(8) :: Amdest(*), Ellhs(*), Flhs(*)
      integer :: ntemp(10000)
!
      integer :: n,i,j,k,l,mi,mj,mk
!
      n = 1
!
      if (NDOFM .gt. 10000) then
         write(NFSOUT,7771)
7771      format(1x,'THE TEMP ARRAY: ntemp MUST ', &
                 'BE INCREASED IN SYMASM')
         IERR = 1
         return
      endif

!      write(*,*) 'in symasm ellhs = '
!      kk=0
!      do ii=1,NDOFM
!        do jj=1,ii
!          kk=kk+1
!          write(*,*) Ellhs(kk)
!        enddo
!        pause
!      enddo

!
! loop thru the dof in the element
!
      do i = 1,NDOFM
!     ------------------
! pull the dof destination vector (for this row)
!
          mi = int(Amdest(i))
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
         do j = 1, i
            ntemp(j) = n + j - 1
         enddo
!
! ALLIANT directives
!vd$ nodepchk
!vd$ nosync
! ARDENT directives
!$doit IVDEP
!wb <
! becuz its symmetric loop thru lower triangular dof for this element
!
         do j = 1, i
!        -------------
! pull the dof destination vector (for this column)
!
            mj = int(Amdest(j))
!wb >
! determine the position in the front (function locls)
!wb   (**note: we may wish to inline this code)
!
!wb             mk = locls(mi,mj)
            k = max0(mi,mj)
            l = min0(mi,mj)
!
            mk = k*(k-1)/2 + l
!wb <
!
! assemble into the front
!wb >
!wb             Flhs(mk) = Flhs(mk) + Ellhs(n)
!wb             n = n + 1
            Flhs(mk) = Flhs(mk) + Ellhs(ntemp(j))
!wb <
!
         enddo
!wb >
         n = n + i
!wb <
!
      enddo
!
!
   end subroutine symasm
