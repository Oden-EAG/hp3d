!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: Primary driver routine for the frontal solver solution modul
! --------------------------------------------------------------
! surfs ===> symmetric / unsymmetric / resolution frontal solver
! --------------------------------------------------------------
! **note: you must have already set up the destination vectors
!          using the prefront routines (sub.surfsp, etc)
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! W   Awrk  : temporary workspace required by surfs
!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
! LATEST REVISION: Feb 2023
!++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
! NAMING CONVENTIONS:
!     AAAAAAAA    Variables in COMMON & PARAMETERS
!     Aaaaaaaa    Variables as ARGUMENTS
!     aaaaaaaa    LOCAL Variables
!         7xxx    FORMAT Statements
!         9xxx    ERROR Handling
!+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
   subroutine surfss (Awrk)
!
      use surfsc1
!
      real(8) :: Awrk(*)
!
! debug print
!
      if (IPFSS1 .eq. 1) then
         write(NFSOUT,7701) ISYM,NUMELM,IRESOL,NRHS,MA,IWRT, &
                            IASSEM,NFSOUT,NICMUL
7701     format(1x,'SURFSS: ISYM,NUMELM,IRESOL,NRHS,MA,IWRT,', &
                           'IASSEM,NFSOUT,NICMUL',/,10x,9i8)
      endif
!
! for no resolution, solve the system of equations
!
      if (IRESOL .eq. 0) call complt (Awrk)
!                        ---------------
!
! with resolution, solve the system of equations
!
      if (IRESOL .eq. 1) call resol (Awrk)
!                        -------------
   end subroutine surfss
