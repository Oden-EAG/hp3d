!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: Assembles the rhs's into the front
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   Amdest:  dof destination vectors
! I   Elrhs :  element rhs
! IO  Frhs  :  front assembled rhs
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
   subroutine semrhs (Amdest, Elrhs, Frhs)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      complex(8) :: Amdest(*),Elrhs(*),Frhs(*)
!
      integer :: n,in,i,j,md
!
      n = 1
!
! loop over the number of rhs's
!
      do 70 in=1,NRHS
!     ---------------
!
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!
! loop thru the dof in the element
!
         do 50 i=1,NDOFM
!        ---------------
! pull the dof destination vector
!
            md = int(Amdest(i))
!
!cwb >
! determine the position in the front (function locr)
!cwb   (**note: we may wish to inline this code)
!
!cwb             j = locr(in,md)
            j = (in-1)*MFW + md
!cwb <
!
            Frhs(j) = Frhs(j) + Elrhs(n)
            n = n + 1
   50    continue
!
   70 continue
!
!
   end subroutine semrhs
