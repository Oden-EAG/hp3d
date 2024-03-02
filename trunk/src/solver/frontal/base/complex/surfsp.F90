!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: Primary driver routine for the frontal solver prefront modul
! --------------------------------------------------------------
! surfs ===> symmetric / unsymmetric / resolution frontal solver
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
!
! I   In : array (length NUMELM) containing the number of nicknames/eleme
!          where the number of nicknames = nodes + multipliers (in gener
!          **note: this must be pre-built by the calling program
!
! I   Ia : array that contains element nick names
!           nicknames are packed as:
!
!            nick = (ndof + NICMUL*nodenum)
!
!            where :   ndof    = number of dof associated w/ this nickna
!                      nodenum = node number or some other unique identi
!                                 reflecting the connectivity or relatio
!                                 to the system of equations
!
!        **note: this must be pre-built by the calling program
!
! O   Ms : minimum storage required for symmetric problem
! O   Mr : minimum stroage required per right hand side
! O   Mu : minimum storage required for unsymmetric problem
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
   subroutine surfsp (In, Ia, Ms, Mu, Mr)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      integer :: In(*), Ia(*)
      integer :: Ms,Mu,Mr
!
      real(8) :: dt,t0,t1
      integer :: i,maxia,ntot
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
! initialize the timers
!
      call zecond(t0)
!
      IERR = 0
      NFN = 0
      MLDEST = 0
!
! determine the maximum number of destination vectors associated w/any 1
!-----------------------------------------------------------------------
! loop thru all the elements
!
      do i = 1,NUMELM
!
! simply pick up the maximum number of nodes+multipliers/element
!
         if (In(i) .gt. MLDEST) MLDEST = In(i)
!
! increment NFN,  NFN = summation of nodes+multipliers/element
!
         NFN = NFN + In(i)
!
     enddo
!
! check if there is enuf space in ia() to perform the prefront
!
      maxia = 2*(NFN+MLDEST)
      if (maxia .gt. MA) then
         write(NFSOUT,7010) maxia,MA
 7010     format(2(/),'ERROR IN SURFSP: NOT ENOUGH ROOM IN ARRAY IA',/, &
               5x,'STORAGE REQUIRED ',i7,/, &
               5x,'STORAGE AVAILABLE ',i7)
         IERR = 1
         return
      endif
!
! calculate the destination vectors
! =================================
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(nfsout,*) 'IN SURFSP NFN,MLDEST = ',NFN,MLDEST
      endif
#endif
      call desvec (In, Ia, Ia(NFN+1), Ia(2*NFN+1), ntot)
!     --------------------------------------------------
!
! set the minimum storage requirements
! ------------------------------------
!
      Mr = IASSEM*MDOF + MFW + 1
      Ms = NUMELM + MLDEST + 2*MDOF + IASSEM*MDOF*(MDOF+1)/2 &
                  + mfw*(MFW+1)/2 + MFW
      Mu = NUMELM + MLDEST + 2*MDOF + IASSEM*MDOF*MDOF + MFW*MFW + MFW
!
! pick up the final time
!

      call zecond(t1)

      dt = t1-t0
!
! print the prefront summary
!----------------------------
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write (NFSOUT,7020) MFW,MDOF,MLDEST
!
        if (ISYM.eq.1 .or. ISYM.eq.4) write(NFSOUT,7030) Ms,Mr
        if (ISYM.eq.2 .or. ISYM.eq.3) write(NFSOUT,7040) Mu,Mr
!
        write(NFSOUT,7050) dt
      endif
!
 7020 format(2(/), 5x,'COMPLETED PREFRONT ',// &
                  10x,'MAXIMUM FRONT WIDTH =  ',i5,/ &
                  10x,'MAXIMUM ELEMENT DOF =  ',i5,/ &
                  10x,'MAXIMUM ELEMENT DESTINATION VECTORS =  ',i2)
 7030 format(10x,'SYMMETRIC PROBLEM MIN. STORAGE =  ',i6,' + ', &
                  i4,'*NRHS '   )
 7040 format(10x,'UNSYMMETRIC PROBLEM MIN. STORAGE =  ',i6, &
                  ' + ',i4,'*NRHS'   )
 7050 format(10x,'TIME IN PREFRONT =   ',f6.3,/)
!
#endif
!
   end subroutine surfsp
