!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! Initiate forward elimination of lhs and rhs followed by backsubstituti
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
!
! W   Awrk    Workspace used by the solver (it is cut and passed below)
!
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
      subroutine complt (Awrk)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      complex(8) :: Awrk(*)
!
      real(8) :: dt,t0,tb,tf
      integer :: n,ia1,iab,iae,iaf,ial,iam,ian,iful,minn
!
      call zecond(t0)
!
      NNEGP = 0
      NPOSP = 0
      IERR = 0
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
! set up the storage cut lengths
!--------------------------------
! MDOF = maximum number of dof for any elem
! n = amount of space required to hold the first few temp arrays cut fro
!      the workspace
!
      n = NUMELM + MLDEST + 2*MDOF
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
         write(NFSOUT,*) 'in complt.f'
         write(NFSOUT,*) 'mdof,numelm,mldest,n',mdof,numelm,mldest,n
      endif
#endif
!
! MELEM = space required to hold the maximum elem lhs + rhs, during asse
! MKF = space required to hold the lhs equations for the max frontwidth
!
! symmetric
      if(ISYM.eq.1 .or. ISYM.eq.4) then
         MELEM = MDOF*(MDOF+1)/2 + MDOF*NRHS
         MKF = MFW*(MFW+1)/2
!
! unsymmetric
      else
         MELEM = MDOF*(MDOF+NRHS)
         MKF = MFW*MFW
      endif
!
! MFWR = space required to hold the lhs & rhs equations for the max fron
!
      MFWR = MKF + MFW*NRHS
!
! no elem assembly space is needed when the data comes from outside
!  (ie: IASSEM=0)
!
      MELEM = IASSEM*MELEM
!
! MW = total for front and individual elem temp storage
!
      MW = MELEM + MFWR
!
! space which is left over is allocated for buffering equations
!
      MBUF = MA - MW - n
      if (MBUF.lt.1000) then
        write(NFSOUT,*) 'COMPLT:TOO LITTLE BUFFER !! MBUF = ',MBUF
        stop 'complt 1'
      endif
!cwb >
! save the value of MBUF & MW to check against during resolution
! *** Note: we cannot let the size of MBUF shrink during resolution
!            because then the stuff that was buffered out will not fit b
!
      if (ISYM .eq. 3) then
         MBUFSV(1,1) = MBUF
         MBUFSV(1,2) = MW
      elseif (ISYM .eq. 4) then
         MBUFSV(2,1) = MBUF
         MBUFSV(2,2) = MW
      endif

!cwb <
!
! open direct access buffer files
! ===============================
! open the LHS storage buffer diskfile
!
      if(IDUMPWR.eq.1)then
        write(NFSOUT,*)'COMPLT: PRETEND OPEN U'
        call zdiriodum('U', 'OPEN', 0, MBUF, Awrk, IERR)
      else
        call zdirio ('U', 'OPEN', 0, MBUF, Awrk, IERR)
!       --------------------------------------------
      endif

!
C open L : The additional eliminated LHS buffer for unsymmetric w/ resol
!
      if (ISYM .eq .3) call zdirio ('L', 'OPEN', 0, MBUF, Awrk, IERR)
!                     -------------------------------------------
      IERR = 10*IERR
      if (IERR .ne. 0) go to 9100
!
! debug print
!
      if (IPRSTR .eq. 1) then
         if(ISYM.eq.1 .or. ISYM.eq.4) then
           write(NFSOUT,7000)
         else
           write(NFSOUT,7010)
         endif
         if(ISYM.eq.2 .or. ISYM.eq.1) write(NFSOUT,7020)
         write(NFSOUT,7030) n, MELEM, MFWR, MBUF, MA
!
 7000    format(2(/), 5x,'SYMMETRIC FORWARD ELIMINATION',/)
 7010    format(2(/), 5x,'UNSYMMETRIC FORWARD ELIMINATION',/)
 7020    format(      5x,'RESOLUTION INACTIVATED',/)
 7030    format(     4x,'            OVER HEAD:',i7,/, &
                     4x,'              ELEMENT:',i7,/, &
                     4x,'                FRONT:',i7,/, &
                     4x,'               BUFFER:',i7,/, &
                     4x,'        TOTAL STORAGE:',i7)
      endif
!
! check that the buffer workspace is large enuf
!   it must be able to at least hold the max frontwidth and NRHS rightha
!
      if(MBUF .lt. MFW+NRHS) then
         IERR = 6
         minn = MA + MFW + NRHS - MBUF
         write(NFSOUT,7060) minn
 7060    format(2(/), 5x,'ERROR IN SURFS: NOT ENOUGH ROOM IN BUFFER', &
                   /,12x,'MINIMUM SIZE OF WORK ARRAY  = ',i7)
         go to 9100
      endif
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
      ial  =  ia1 + NUMELM
!
! set pointer to temp storage for dof destination vectors         (amdes
!
      iam  =  ial + MLDEST
!
! set pointer to temp storage for the transfer loc for active dof
!   in the front after eliminating                                (andes
!
      ian  =  iam + MDOF
!
! set pointer to temp storage for elem lhs + rhs, during assembly   (ele
!
      iae  =  ian + MDOF
!
! set pointer to temp storage for lhs & rhs equations for front     (frn
!
      iaf  =  iae + MELEM
!
! set pointer to temp storage for equation buffering,
!   this is all remaining workspace                                  (bu
!
      iab  =  iaf + MFWR
!
! debug print
!
      if (IPFSST .eq. 1) then
         write(NFSOUT,7701) NUMELM,MLDEST,MDOF,NRHS,MELEM,MKF,MFWR
 7701     format(1x,'COMPLT: NUMELM,MLDEST,MDOF,NRHS,MELEM,MKF,MFWR', &
                 /,10x,8i10)
         write(NFSOUT,7702) MBUF,MA,MW,ia1,ial,iam,ian,iae,iaf,iab
 7702     format(1x,'COMPLT: MBUF,MA,MW,ia1,ial,iam,ian,iae,iaf,iab', &
                 /,10x,10i10)
      endif

!
! forward elimination of all elements  (lhs and rhs)
! ====================================
!
      call frwcp (Awrk(ia1), Awrk(ial), Awrk(iam), Awrk(ian), &
                  Awrk(iae), Awrk(iaf), Awrk(iab))
!     --------------------------------------------------------
!
      call zecond(tf)
      dt  =  tf-t0
!
      if(IPRSTR .eq. 1)then
         write(NFSOUT,7080) dt
         iful = IFU*IDUMP + IFL
         if(iful .ne. 0) write(NFSOUT,7100) iful
 7080     format(10x,'TIME IN FORWARD ELIMINATION:',f9.3)
 7100     format(10x,'BUFFER DUMPS:',i4)
      endif
!
      if (IERR .ne. 0) go to 9100
!
      if (NRHS .ne. 0) then

!
!
! backward substitution of all elements
! ====================================
!
         call bckwrd (Awrk(ia1), Awrk(ial), Awrk(iam), Awrk(ian), &
                      Awrk(iae), Awrk(iaf), Awrk(iab), Awrk(iab))
!        --------------------------------------------------------
!
         call zecond(tb)
         dt  =  tb - tf
         if(IPRSTR .eq. 1) write (NFSOUT,7120) dt
 7120     format(10x,'TIME IN BACK SUBSTITUTION:',f9.5)
      endif
!
! normal exit
!------------
!
      return
!
! error exit
! ==========
!
 9100 continue
!
!
      end subroutine complt
