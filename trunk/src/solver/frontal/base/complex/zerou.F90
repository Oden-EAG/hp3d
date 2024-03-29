!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: If necessary, zero out for new unsymmetric lhs equations
!            in the front
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
!
! IO  Flhs  : the space which holds the lhs equations in the front
!          x  note: for unsymmetric, lhs is stored as:
!               a11 a12 a13 ...
!               a21 a22 a23 ...  ==>  [a11,a12,a13,..., a21,a22,a23,...]
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
   subroutine zerou (Flhs)
!
      use surfsc2
!
      implicit none
!
      complex(8) :: Flhs(*)
!
      integer :: i,j,mi,mj,mk
!
      complex(8), parameter :: zero = (0.d0, 0.d0)
!
!
! if the new frontwidth is equal to the old frontwidth
!   then there is nothing to zero
!
      if(NFW .eq. LFW) return
!
      mi = LFW*NFW + 1
      mj = NFW*NFW
      mk = LFW*LFW + 1
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
      do i = mi,mj
        Flhs(i) = zero
      enddo
!
      if (LFW .eq. 0) return
!
      mj = NFW - LFW
!
      do i = 1,LFW
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
         do j = 1,mj
            mi = mi - 1
            Flhs(mi) = zero
         enddo
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
! also shove to the end the RHS elimination terms
!
         do j = 1,LFW
            mi = mi - 1
            mk = mk - 1
            Flhs(mi) = Flhs(mk)
         enddo
!
      enddo
!
!
   end subroutine zerou
