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
      complex(8) :: Frhs(*)
!
      integer :: i,ia,in,m
!
      complex(8), parameter :: zero = (0.d0, 0.d0)
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
! if the new frontwidth is equal to the old frontwidth
!   then there is nothing to zero
!
      if (LFW .eq. NFW) return
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(nfsout,*) 'ZEROR: NRHS,MFW,NFW,LFW = ',NRHS,MFW,NFW,LFW
      endif
#endif
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
            Frhs(ia+i) = zero
         enddo
!
      enddo
!
!
   end subroutine zeror
