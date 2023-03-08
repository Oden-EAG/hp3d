      subroutine surfsp (In, Ia, Ms, Mu, Mr)
c
#include "syscomz.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: Primary driver routine for the frontal solver prefront modul
c --------------------------------------------------------------
c surfs ===> symmetric / unsymmetric / resolution frontal solver
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c
c I   In : array (length numel) containing the number of nicknames/eleme
c          where the number of nicknames = nodes + multipliers (in gener
c          **note: this must be pre-built by the calling program
c
c I   Ia : array that contains element nick names
c           nicknames are packed as:
c
c            nick = (ndof + NICMUL*nodenum)
c
c            where :   ndof    = number of dof associated w/ this nickna
c                      nodenum = node number or some other unique identi
c                                 reflecting the connectivity or relatio
c                                 to the system of equations
c
c        **note: this must be pre-built by the calling program
c
c O   Ms : minimum storage required for symmetric problem
c O   Mr : minimum stroage required per right hand side
c O   Mu : minimum storage required for unsymmetric problem
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
#include "surfsc1.blk"
#include "surfsc2.blk"
c
      integer :: In(*), Ia(*)
      real(8) :: t0,t1
c
#if DEBUG_MODE
      integer :: iprint=0
#endif
c
c initialize the timers
c
      call zecond(t0)
c
      IERR = 0
      NFN = 0
      MLDEST = 0
c
c determine the maximum number of destination vectors associated w/any 1
c-----------------------------------------------------------------------
c loop thru all the elements
c
      do 10 i = 1,NUMELM
c
c simply pick up the maximum number of nodes+multipliers/element
c
         if (In(i) .gt. MLDEST) MLDEST = In(i)
c
c increment NFN,  NFN = summation of nodes+multipliers/element
c
         NFN = NFN + In(i)
c
   10 continue
c
c check if there is enuf space in ia() to perform the prefront
c
      maxia = 2*(NFN+MLDEST)
      if (maxia .gt. MA) then
         write(NFSOUT,7010) maxia,MA
 7010     format(2(/),'ERROR IN SURFSP: NOT ENOUGH ROOM IN ARRAY IA',/,
     .         5x,'STORAGE REQUIRED ',i7,/,
     .         5x,'STORAGE AVAILABLE ',i7)
         IERR = 1
         return
      endif
c
c calculate the destination vectors
c =================================
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(nfsout,*) 'IN SURFSP NFN,MLDEST = ',NFN,MLDEST
      endif
#endif
      call desvec (In, Ia, Ia(NFN+1), Ia(2*NFN+1), ntot)
c     --------------------------------------------------
c
c set the minimum storage requirements
c ------------------------------------
c
      Mr = IASSEM*MDOF + MFW + 1
      Ms = NUMELM + MLDEST + 2*MDOF + IASSEM*MDOF*(MDOF+1)/2
     .            + mfw*(MFW+1)/2 + MFW
      Mu = NUMELM + MLDEST + 2*MDOF + IASSEM*MDOF*MDOF + MFW*MFW + MFW
c
c pick up the final time
c

      call zecond(t1)

      dt = t1-t0
c
c print the prefront summary
c----------------------------
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write (NFSOUT,7020) MFW,MDOF,MLDEST
c
        if (ISYM.eq.1 .or. ISYM.eq.4) write(NFSOUT,7030) Ms,Mr
        if (ISYM.eq.2 .or. ISYM.eq.3) write(NFSOUT,7040) Mu,Mr
c
        write(NFSOUT,7050) dt
      endif
c
 7020 format(2(/), 5x,'COMPLETED PREFRONT ',//
     .            10x,'MAXIMUM FRONT WIDTH =  ',i5,/
     .            10x,'MAXIMUM ELEMENT DOF =  ',i5,/
     .            10x,'MAXIMUM ELEMENT DESTINATION VECTORS =  ',i2)
 7030 format(10x,'SYMMETRIC PROBLEM MIN. STORAGE =  ',i6,' + ',
     .            i4,'*NRHS '   )
 7040 format(10x,'UNSYMMETRIC PROBLEM MIN. STORAGE =  ',i6,
     .            ' + ',i4,'*NRHS'   )
 7050 format(10x,'TIME IN PREFRONT =   ',f6.3,/)
c
#endif
c
      end subroutine surfsp
