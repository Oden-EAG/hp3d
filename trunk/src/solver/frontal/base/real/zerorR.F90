!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: If necessary, zero out for new rhs equations in the front
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! IO  Frhs  : the space which holds the rhs equations in the front
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
   subroutine zeror (Frhs)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      real(8) :: Frhs(*)
!
      integer :: i,ia,in,m
!
      real(8), parameter :: dzero = 0.0d0
!
!  ...if the new frontwidth is equal to the old frontwidth
!     then there is nothing to zero
!
      if (LFW .eq. NFW) return
!
      do in=1,NRHS
!
         ia = (in-1) * MFW
         m = LFW + 1
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
         do i = m,NFW
            Frhs(ia+i) = dzero
         enddo
!
      enddo
!
!
   end subroutine zeror
