!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! Forward elimination of both lhs and rhs
!  calls solin1 for dest. vectors
!  calls solin2 for lhs and rhs's
!
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! W   Alelm   length of the current front as we process each Elem
! W   Aldest  nodal destination vectors
! W   Amdest  dof destination vectors
! W   Andest  destination in the front where the remaining values
!              in the front will transfer to as each exiting dof is elim
! W   Elem    Element lhs + rhs
! W   Frnt    equations in the front
! W   Buf     eliminated equations
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
   subroutine frwcp (Alelm, Aldest, Amdest, Andest, Elem, Frnt, Buf)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      complex(8) :: Aldest(*),Amdest(*),Andest(*),Elem(*),Frnt(*), &
                    Buf(*),Alelm(*)
!
      integer :: i,id,ie,iel,ifg,iunsr,j,jerr,lenf
      integer :: m,md,mke,n,ne,negiel,newl,newu,numdes
!
!  ...test......test......test......test......test......test......test...
!!!      integer :: ndest1(100)
!  ...test......test......test......test......test......test......test...
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!
      IFU = 0
      IFL = 0
      NRHSF = NRHS
      IU = 1
!
      IL = MBUF
!
      NFW = 0
      LFW = 0
!
      IDUMP = 0
      LENU  = 0
!
      iunsr = 0
      if (ISYM.eq.3) iunsr=1
!
! loop over the Elements
! **********************
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(nfsout,*) 'FRWCP: NUMELM = ',NUMELM
      endif
#endif
      do iel=1,NUMELM
!     ********************


!
         negiel = -iel
!
! pull up the Elements nodal destination vectors (Aldest)
!
         call solin1 (iel, numdes, Aldest)
#if HP3D_DEBUG
         if (iprint.eq.1) write(*,*) 'iel,numdes',iel,numdes
#endif
!        ------------------------------


!  ...test......test......test......test......test......test......test...
!  ...test......test......test......test......test......test......test...
!         read(30,*)iel1
!         read(30,*)numdes1
!         read(30,*)(ndest1(i),i=1,numdes1)
!
!         if(numdes.ne.numdes1)then
!           write(*,*)'FRWCP: iel,numdes,numdes1=',iel,numdes,numdes1
!           stop
!         endif
!         do i=1,numdes
!            if(ndest1(i) .ne. nint(dreal(aldest(i))))then
!              write(*,*)'FRWCP:IEL,I=',iel,i
!              write(*,*)ndest1(i),nint(dreal(aldest(i)))
!              stop
!            endif
!            enddo
!  ...test......test......test......test......test......test......test...
!  ...test......test......test......test......test......test......test...


         if(numdes .eq. 0) goto 200
!
! convert the nodal destination vectors to dof destination vectors (Amde
!  also determine the number of dof to eliminate from the front (ne)
!  and the transfer loc for active dof in the front after eliminating (A
!
         call dest (numdes, Aldest, ne, Amdest, Andest)
!        -----------------------------------------


!
! check that the frontwidth has not reduced, if so, dont let it
!  this is necessary especially for indexing for unsymmetric k
!
         if (LFW .gt. NFW) NFW = LFW
         KFW = NFW
         mke = MKF
!
! zero out for any new dof in the front
!
! symmetric
!
         if (ISYM.eq.1 .or. ISYM.eq.4) then
            if(IASSEM .eq. 1) mke = NDOFM*(NDOFM+1)/2
            call zeros (Frnt)
!           -----------------
! unsymmetric
!
         else
            if(IASSEM .eq. 1) mke = NDOFM**2
            call zerou (Frnt)
!           -----------------
         endif
!
         if (NRHS .ne. 0) call zeror (Frnt(MKF+1))
!                         ------------------------
!
! build the Element rhs & lhs
!----------------------------
! pointers into Elem(): (1: points to lhs  ;  mke+1: points to rhs)
! flag for solin2 (ifg=1: get rhs ; =2: get lhs as well)
!

         ifg = 2
         md = int(Amdest(1))
         call solin2 (iel, ifg, NDOFM, NRHS, md, Elem, Elem(mke+1))
!        --------------------------------------------------------------


! debug print
!
         if (IPFSST .eq. 1 .or. IPFSST.eq.negiel) then
            write (NFSOUT,7701) iel,NRHS,NDOFM,NFW,KFW,LFW,NE,mke,mkf
 7701        format(1x,'FRWCP: iel,NRHS,NDOFM,NFW,KFW,LFW,NE,mke,mkf', &
                    /,10x,10i10)
         endif
         if (IPFSLH .eq. 1 .or. IPFSLH.eq.negiel) then
            write (NFSOUT,7702) (Elem(i),i=1,mke)
 7702      format(1x,'FRWCP: Elem: Element LHS',/,(10x,1p,10e12.5))
         endif
         if (IPFSRH .eq. 1 .or. IPFSRH.eq.negiel) then
            write (NFSOUT,7703) ((Elem(mke+i),i=1,NDOFM),j=1,NRHS)
 7703      format(1x,'FRWCP: Elem: Element RHS',/,(10x,1p,10e12.5))
         endif
!
! assemble the Element rhs & lhs into the front
! ---------------------------------------------

         if (IASSEM .ne. 0) then
!wb >
! to save having to do abs() repeatedly, do it now
!  **note: sub.dest is recalled in sub.bckwrd so this wont screw it up
!          BUT the negative amdest info is now lost herein
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!


            do id = 1,NDOFM
!wr10.07.99
!!!!!          Amdest(id) = Cdabs(Amdest(id))
               Amdest(id) =   abs(Amdest(id))
            enddo

!wb <
!
! symmetric
!

            if (ISYM.eq.1 .or. ISYM.eq.4) then
               call symasm (Amdest, Elem, Frnt)
!              --------------------------------
! unsymmetric
!
            else

               call unsasm (Amdest, Elem, Frnt)
!              --------------------------------
            endif
!wb >
            if (IERR .ne. 0) return
!wb <
!
            if(NRHS .ne. 0) call semrhs (Amdest, Elem(mke+1), &
                                                 Frnt(MKF+1))
!                           ---------------------------------
! debug print
!
            if (IPFSLF .eq. 1 .or. IPFSLF.eq.negiel) then
               lenf = (NFW*(NFW+1))/2
               if (ISYM.eq.2 .or. ISYM.eq.3) lenf = NFW*NFW
               write (NFSOUT,7704) (Frnt(i),i=1,lenf)
 7704          format(1x,'FRWCP: Frnt: Assembled LHS',/, &
                       (10x,1p,5(2e12.5,2x)))
               call pause
            endif
            if (IPFSRF .eq. 1 .or. IPFSRF.eq.negiel) then
               write (NFSOUT,7705) ((Frnt(MKF+i),i=1,NFW),j=1,NRHS)
 7705          format(1x,'FRWCP: Frnt: Assembled RHS',/, &
                      (10x,1p,5(2e12.5,2x)))
               call pause
            endif
!
         endif
!
! eliminate those equations possible from the current front
! **********************************************************

            if (IPFSRF .eq. 1 .or. IPFSRF.eq.negiel) then
              write(*,*) 'frwcp: ne = ',ne
              call pause
            endif

!
         do ie = 1,ne
!        ****************
!
! check if the equation Buffers are going to get full
!   and dump if necessary
! -------------------------------------------------
! *note: Ubuf runs forward into the workspace Buff(1-->MBUF)
!        Lbuf runs backward into the workspace Buff(1<--MBUF)
!
!
            newu = IU + NFW + NRHS - 1
!wb >
!wb             newl = IL + 1 - NFW * ISYM/3
            newl = IL + 1 - NFW * iunsr
!wb <



!
            if (newu .ge. newl) then
!
! increment the record counter
!
               IFU = IFU + 1
!
! dump the Buffer
!
               if(IDUMPWR.eq.1)then
                 write(nfsout,*) &
                 'FRWCP: TOO LITTLE SPACE TO AVOID DUMPS'
                 stop 'frwcp 1'
               endif
               call zdirio('U', 'WRITE', IFU, IU-1, Buf, jerr)
!              -----------------------------------------------
               IDUMP = 1
!
! error out
!
               IERR = 10*jerr
               if(jerr .ne. 0) return
!
! set the pointer to the beginning of the equation Buffer
!
               iu = 1
!
! for unsymmetric w/ resolution
!
               if (ISYM .eq. 3) then
                  IFL = IFL + 1
                  call zdirio ('L','WRITE',IFL,MBUF-IL,Buf(IL+1),jerr)
!                 ---------------------------------------------------
                  IL = MBUF
                  IERR = 10*jerr
                  if (jerr .ne. 0) return
               endif
!
! debug print
!
               if (IPFSZD .eq. 1) then
                  write (NFSOUT,7706)  IFU,IFL,IL,IU,MBUF
 7706              format(1x,'FRWCP: IFU,IFL,IL,IU,MBUF',5i8)
               endif
!
! debug print
!
               if (IPFSBF .eq. 1) then
                  write (NFSOUT,7711) (Buf(i),i=1,40)
 7711              format(1x,'FRWCP: First 40 items in Buf:',/, &
                           (10x,1p,10e12.5))
!
                  if (ISYM .eq. 3) write (NFSOUT,7712) &
                                   (Buf(i),i=IL+1,il+40)
 7712              format(1x,'FRWCP: First 40 items in Buf(IL+1):',/, &
                          (10x,1p,10e12.5))
!
               endif
!
            endif
!
! eliminate this equation lhs from the current front
!  ----------------------------------------------
            m = IU
!
! symmetric
!
            if (ISYM.eq.1 .or. ISYM.eq.4) then
               call symelm (iel, Andest(ie), Frnt, Buf(IU))
!              --------------------------------------------
! unsymmetric
!
            else
               call unselm (iel, Andest(ie), Frnt, Buf(IU))
!              --------------------------------------------
            endif
!
! debug print
!
            if (IPFSLE .eq. 1 .or. IPFSLE.eq.negiel) then
              lenf = (NFW*(NFW+1))/2
              if (ISYM.eq.2 .or. ISYM.eq.3) lenf = NFW*NFW
              write (NFSOUT,7707) iel,Andest(ie),(Frnt(i),i=1,lenf)
 7707         format(1x,'FRWCP: iel,Andest(ie):',i8,1p,(2e12.4,2x), &
                     /,10x,'Frnt: Eliminated LHS',/, &
                     (10x,1p,5(2e12.5,2x)))
              call pause
            endif
!
! increment the equation Buffer pointer
! ------------------------------------
!
            IU = IU + NRHS + NFW
!
            if (IERR .eq. 2) then
               write (NFSOUT,7000) iel
 7000           format(2(/), 5x,'ERROR IN FRWCP:', &
                      ' ZERO PIVOT IN ELEMENT:',i4)
               return
            endif
!
! eliminate this equation rhs from the current front
!  ----------------------------------------------
            if(NRHS  .ne.  0) then
!
! symmetric
!
               if (ISYM.eq.1 .or. ISYM.eq.4) then
                  call elmrhs (Andest(ie),   1, Frnt(MKF+1), Buf(m), &
                               Buf(m+NFW))
!                 ---------------------------------------------------
! unsymmetric
!
               else
                  call elmrhs (Andest(ie), KFW, Frnt(MKF+1), Frnt(NFW), &
                               Buf(m+NFW))
!                 -----------------------------------------------------
               endif
!
! debug print
!
               if (IPFSRE .eq. 1 .or. IPFSRE.eq.negiel) then
                  write (NFSOUT,7708) iel,Andest(ie), &
                                      ((Frnt(MKF+i),i=1,NFW),j=1,NRHS)
 7708              format(1x,'FRWCP:  iel,Andest(ie):',i8,1p, &
                          (2e12.4,2x),/,10x,'Frnt: Eliminated RHS',/, &
                          (10x,1p,5(2e12.5,2x)))
               endif
!
            endif
!
! for unsymmetric w/ resolution
!   we must copy the the k(i,m)/k(m,m) terms into the Lbuf Buffer space
!   it may be dumped to diskfile for resolution
!   (see sub.unselm and sub.elmrhs for further explanation)
!   *note: Ubuf runs forward into the workspace Buff(1-->MBUF)
!          Lbuf runs backward into the workspace Buff(1<--MBUF)
!
            if (ISYM .eq. 3) then
               m = NFW
               n = NFW-1
               n = max0(n,1)
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
!
               do j = 1,n
                  Buf(IL) = Frnt(m)
                  IL = IL - 1
                  m = m + KFW
               enddo
!
            endif
!
! reduce the front by the equation we just reduced
!
            NFW = NFW - 1
!           -------------
         enddo
! ***************
! store the final length of the current front (LFW) into Alelm
!
#if HP3D_DEBUG
         if (iprint.eq.1) then
           write(*,*) 'frwcp: NFW = ',NFW
         endif
#endif
         LFW = NFW
         Alelm(iel) = LFW
!
! for unsymmetric w/ equations being eliminated
!   we must copy the the k(i,m)/k(m,m) terms forward
!   (see sub.unselm and sub.elmrhs for further explanation)
!
!
         if ((ISYM.eq.2 .or. ISYM.eq.3) .and. ne.ne.0) then
!
! n = the beginning frontwith
! m = the final frontwidth
!
            n = KFW
            m = NFW + 1
!
            do i = 2,NFW
!
! ALLIANT directives
!vd$ select (vector)
! ARDENT directives
!$doit VBEST
!
!
               do j = 1,NFW
                  Frnt(m) = Frnt(n+j)
                  m = m + 1
               enddo
!
               n = n + KFW
            enddo
!
         endif
#if HP3D_DEBUG
         if (iprint.eq.1) then
           write(*,*) 'frwcp: at the end of the loop'
         endif
#endif
!
  200 enddo
! ************
!
! dump the final Buffers to diskfiles
!  (or the only Buffers, full or not, for resolution)
!---------------------------------------------------------------
!
      IFU = IFU + 1
!wb >
!wb       call zdirio ( 'u', 'write', IFU, iu-1, Buf, jerr)
!wb !     --------------------------------------------------
!wb
      jerr = 0
!
! if we are not going to do a resolution (ISYM=3 OR 4)
!   and we have not had to dump any Buffers,  then skip this write
!


      if (ISYM.eq.4 .or. IFU.gt.1 .or. ISYM.eq.3) then

         IDUMP = 1

         if(IFU.le.1 .and. ISYM.ne.3) &
           then
!           ...ISYM=4 is the exclusive reason to dump:
               if(IDUMPWR.eq.1) &
                 then
                     write(nfsout,*) &
                     'FRWCP: PRETEND WRITING, if ISYM=4 !!!'
                     call zdiriodum('U', 'WRITE', IFU, IU-1, Buf, jerr)
                 else
                     call zdirio   ('U', 'WRITE', IFU, IU-1, Buf, jerr)
                 endif
           else
               if(IDUMPWR.ne.1) &
                 then
                     call zdirio   ('U', 'WRITE', IFU, IU-1, Buf, jerr)
!                    -------------------------------------------------
                 else
                     write(nfsout,*)'FRWCP: AVOIDING DUMP IMPOSSIBLE'
                     stop 'frwcp 2'
                 endif
           endif
!
      else
         LENU = IU - 1
      endif
!wb <
      IERR = 10*jerr
      if(jerr .ne. 0) return
!
      if (ISYM .eq. 3) then
         IFL = IFL+1
         call zdirio ('L', 'WRITE', IFL, MBUF-IL, Buf(IL+1), jerr)
!        ----------------------------------------------------------
         IERR=10*jerr
         if(jerr .ne. 0) return
      endif
!
! debug print
!
      if (IPFSBF .eq. 1) then
         write (NFSOUT,7709) (Buf(i),i=1,40)
 7709     format(1x,'FRWCP: First 40 items in Buf: ',/, &
                           (10x,1p,10e12.5))
!
         if (ISYM .eq. 3) write (NFSOUT,7710) (Buf(i),i=il+1,il+40)
 7710                   format(1x,'FRWCP: First 40 items in Buf(il+1):', &
                               /,(10x,1p,10e12.5))
!
      endif
!
!
   end subroutine frwcp
