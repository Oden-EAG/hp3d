c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: If necessary, zero out for new rhs equations in the front
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c IO  Frhs  : the space which holds the rhs equations in the front
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
      subroutine zeror (Frhs)
c
      use surfsc1
      use surfsc2
c
      implicit none
c
      complex(8) :: Frhs(*)
c
      integer :: i,ia,in,m
c
      complex(8), parameter :: zero = (0.d0, 0.d0)
c
#if DEBUG_MODE
      integer :: iprint
      iprint=0
#endif
c
c if the new frontwidth is equal to the old frontwidth
c   then there is nothing to zero
c
      if (LFW .eq. NFW) return
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(nfsout,*) 'ZEROR: NRHS,MFW,NFW,LFW = ',NRHS,MFW,NFW,LFW
      endif
#endif
      do 120 in  =1,NRHS
c
         ia = (in-1) * MFW
         m = LFW + 1
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
         do 20 i = m,NFW
            Frhs(ia+i) = zero
20       continue
c
120   continue
c
c
      end subroutine zeror
