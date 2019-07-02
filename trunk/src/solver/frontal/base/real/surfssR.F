      subroutine surfss (Awrk)
c
#include "syscom.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: Primary driver routine for the frontal solver solution modul
c --------------------------------------------------------------
c surfs ===> symmetric / unsymmetric / resolution frontal solver
c --------------------------------------------------------------
c **note: you must have already set up the destination vectors
c          using the prefront routines (sub.surfsp, etc)
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c W   Awrk  : temporary workspace required by surfs
c*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
c AUTHOR:                             LATEST REVISION:
c++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
c NAMING CONVENTIONS:
c     AAAAAAAA    Variables in COMMON & PARAMETERS
c     Aaaaaaaa    Variables as ARGUMENTS
c     aaaaaaaa    LOCAL Variables
c         7xxx    FORMAT Statements
c         9xxx    ERROR Handling
c+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
c
#include "surfsc1.blk"
c
      dimension Awrk(*)
c
c debug print
c
      if (IPFSS1 .eq. 1) then
         write(NFSOUT,7701) ISYM,NUMELM,IRESOL,NRHS,MA,IWRT,
     .                       IASSEM,NFSOUT,NICMUL
7701      format(1x,'SURFSS:  ISYM,NUMELM,IRESOL,NRHS,MA,IWRT,',
     .           'IASSEM,NFSOUT,NICMUL',/,10x,9i8)
      endif
c
c for no resolution, solve the system of equations
c
      if (IRESOL .eq. 0) call complt (Awrk)
c                        ---------------
c
c with resolution, solve the system of equations
c
      if (IRESOL .eq. 1) call resol (Awrk)
c                        -------------
      return
      end
