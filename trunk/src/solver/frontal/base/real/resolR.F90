!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! initiate forward elimination of rhs only,  followed by backsubstitutio
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! W   Awrk  : workspace used by the solver (it is cut and passed below)
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
   subroutine resol (Awrk)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      real(8) :: Awrk(*)
!
      real(8) :: dt,t0,tb,tf
      integer :: ia1,iabf,iabr,iae,iaf,ial,iam,ian
      integer :: machk,minn,n
!
      call zecond(t0)
!
      NNEGP = -1
      NPOSP = -1
!
      IFU = -1
      IFL = -1
!
      IERR = 0
      IFB = 0
!
! check for invalid circumstances
! -------------------------------
      if (ISYM.eq.2 .or. ISYM.eq.1) then
         IERR=3
         write(NFSOUT,7020)
 7020    format(2(/),5x, &
           'ERROR IN SURFS: ATTEMPTED RESOLUTION WITH', &
           ' RESOLUTION INACTIVATED (IE: ISYM=2 OR 1)')
         return
      endif
!
      if (NRHS .eq. 0) then
         IERR=4
         write (NFSOUT,7040)
 7040    format(2(/), 5x, &
           'ERROR IN SURFS: ATTEMPTED RESOLUTION WITH 0 RHS')
         return
      endif
!
! n = amount of space required to hold the first few temp arrays cut fro
!      the workspace
!
      n = NUMELM + MLDEST + 2*MDOF
!
!
!wb >
!  Pull up the old values of MBUF & MW to check against for resolution
! *** Note: we cannot let the size of MBUF shrink during resolution
!            because then the stuff that was buffered out will not fit b
!
      if (ISYM .eq. 3) then
         MBUF = MBUFSV(1,1)
         MW = MBUFSV(1,2)
      elseif (ISYM .eq. 4) then
         MBUF = MBUFSV(2,1)
         MW = MBUFSV(2,2)
      endif
!
! check against the current size of the workspace
!
      machk = MBUF + MW + n
      if (MA .lt. machk) then
         IERR=6
         write(NFSOUT,7027)
 7027    format(2(/),5x, &
           'ERROR IN SURFS: ATTEMPTED RESOLUTION WITH', &
           ' LESS WORKSPACE THAN WAS PREVIOUSLY ALLOCATED')
         return
      endif
!
!wb <
!
!
! set up storage cut lengths
!---------------------------
! assembly space cut length for elem rhs assembly
!
      MELEM = MDOF*NRHS*IASSEM
!
! cut length for the eleminated rhs equations
!
      MFWR = MFW*NRHS
!
! space which is available for buffering rhs equations
!
      MB = MW - MELEM - MFWR
!
! open direct access buffer files
!
      call zdirio ('B', 'OPEN', 0 , MB, Awrk, IERR)
!     -----------------------------------------
! error out
!
      IERR = 10*IERR
      if (IERR .ne. 0) goto 9100
!
! debug print
!
      if(IPRSTR.eq.1) then
         if(ISYM .eq. 4) then
!           write(NFSOUT,7060) NRHS
         else
!           write(NFSOUT,7080) NRHS
         endif
!        write(NFSOUT,7100) n,MELEM,MFWR,MB,MBUF,MA
!
 7060    format(2(/), 5x,'SYMMETRIC RESOLUTION WITH',i3,' RHS',/)
 7080    format(2(/), 5x,'UNSYMMETRIC RESOLUTION WITH',i3,' RHS',/)
 7100    format(     4x,'          OVER HEAD:',i7,/, &
                     4x,'            ELEMENT:',i7,/, &
                     4x,'              FRONT:',i7,/, &
                     4x,'         RHS BUFFER:',i7,/, &
                     4x,'         LHS BUFFER:',i7,/, &
                     4x,'      TOTAL STORAGE:',i7)
      endif
!
! not enuf space to buffer even 1 equation
!
      if (MB .lt. NRHS) then
         IERR = 5
         minn = max0(0,MB)
         write (NFSOUT,7120) minn
 7120    format(2(/),5x,'ERROR: TOO MANY RHS', &
                   /,12x,'MAXIMUM NUMBER OF RHS =',i2)
         goto 9100
      endif
!
!
! set pointers into the workspace (note: relative pointers)
! ===============================
! set pointer for temp storage of the length of the current front
!   for each elem                                                  (alel
!
      ia1 = 1
!
! set pointer to temp storage for nodal destination vectors       (aldes
!
      ial = ia1 + NUMELM
!
! set pointer to temp storage for dof destination vectors         (amdes
!
      iam = ial + MLDEST
!
! set pointer to temp storage for the transfer loc for active dof
!   in the front after eliminating                                (andes
!
      ian = iam + MDOF
!
! set pointer to temp storage for elem lhs + rhs, during assembly   (ele
!
      iae = ian + MDOF
!
! set pointer to temp storage for rhs equations for front           (frn
!
      iaf = iae + MELEM
!
! set pointer to temp storage for buffered rhs equations            (bbu
!
      iabr = iaf + MFWR
!
! set pointer to temp storage for equation buffering,
!   this is all remaining workspace                                  (bu
!
      iabf = iabr + MB
!
!
! forward elimination of all elements  (rhs only)
! ====================================
!
      call frwrs (Awrk(ia1), Awrk(ial), Awrk(iam), Awrk(ian), &
                  Awrk(iae), Awrk(iaf), Awrk(iabr),Awrk(iabf))
!     ---------------------------------------------
!
      call zecond(tf)
      dt = tf-t0
!
      if(IPRSTR .eq. 1) then
         write(NFSOUT,7140) dt
         if(IFB .gt. 0) write(NFSOUT,7160) IFB
 7140     format(10x,'TIME IN FORWARD ELIMINATION:',f9.3)
 7160     format(10x,'RHS BUFFER DUMPS:',i4)
      endif
!
      if (IERR .ne. 0) goto 9100
!
! backward substitution of all elements
! ====================================
!
      call bckwrd ( Awrk(ia1), Awrk(ial), Awrk(iam), Awrk(ian), &
                    Awrk(iae) ,Awrk(iaf), Awrk(iabr),Awrk(iabf))
!      ----------------------------------------------
!
      call zecond(tb)
      dt = tb - t0
      if(IPRSTR .eq. 1) write(NFSOUT,7180) dt
 7180   format(10x,'TIME IN BACK SUBSTITUTION:',f9.3)
!
! normal exit
!------------
!
      return
!
! error exit
! ----------
 9100 continue
!
!
   end subroutine resol
