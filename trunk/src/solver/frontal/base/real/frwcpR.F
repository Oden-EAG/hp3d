c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c Forward elimination of both lhs and rhs
c  calls solin1 for dest. vectors
c  calls solin2 for lhs and rhs's
c
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c W   Alelm   length of the current front as we process each Elem
c W   Aldest  nodal destination vectors
c W   Amdest  dof destination vectors
c W   Andest  destination in the front where the remaining values
c              in the front will transfer to as each exiting dof is elim
c W   Elem    Element lhs + rhs
c W   Frnt    equations in the front
c W   Buf     eliminated equations
c
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
      subroutine frwcp (Alelm, Aldest, Amdest, Andest, Elem, Frnt, Buf)
c
      use surfsc1
      use surfsc2
c
      implicit none
c
      real(8) :: Alelm(*), Aldest(*), Amdest(*), Andest(*)
      real(8) :: Elem(*),Frnt(*),Buf(*)
c
      integer :: i,id,ie,iel,ifg,iunsr,j,jerr
      integer :: lenf,m,mke,n,ne,negiel,newl,newu,numdes
c
      IFU = 0
      IFL = 0
      NRHSF = NRHS
      IU = 1
c
      IL = MBUF
c
      NFW = 0
      LFW = 0
c
      IDUMP = 0
      LENU  = 0
c
      iunsr = 0
      if (ISYM.eq.3) iunsr=1
c
c loop over the Elements
c **********************
      do 200 iel=1,NUMELM
c     ********************
c
         negiel = -iel
c
c pull up the Elements nodal destination vectors (Aldest)
c
         call solin1 (iel, numdes, Aldest)
c        ------------------------------
         if(numdes .eq. 0) go to 200
c
c convert the nodal destination vectors to dof destination vectors (Amde
c  also determine the number of dof to eliminate from the front (ne)
c  and the transfer loc for active dof in the front after eliminating (A
c
         call dest (numdes, Aldest, ne, Amdest, Andest)
c        -----------------------------------------
c
c check that the frontwidth has not reduced, if so, dont let it
c  this is necessary especially for indexing for unsymmetric k
c
         if (LFW .gt. NFW) NFW = LFW
         KFW = NFW
         mke = MKF
c
c zero out for any new dof in the front
c
c symmetric
c
         if (ISYM.eq.1 .or. ISYM.eq.4) then
            if(IASSEM .eq. 1) mke = NDOFM*(NDOFM+1)/2
            call zeros (Frnt)
c           -----------------
c unsymmetric
c
         else
            if(IASSEM .eq. 1) mke = NDOFM**2
            call zerou (Frnt)
c           -----------------
         endif
c
         if (NRHS .ne. 0) call zeror (Frnt(MKF+1))
c                         ------------------------
c
c build the Element rhs & lhs
c----------------------------
c pointers into Elem(): (1: points to lhs  ;  mke+1: points to rhs)
c flag for solin2 (ifg=1: get rhs ; =2: get lhs as well)
c
         ifg = 2
         call solin2 (iel, ifg, NDOFM, NRHS, Amdest, Elem, Elem(mke+1))
c        --------------------------------------------------------------
c debug print
c
         if (IPFSST .eq. 1 .or. IPFSST.eq.negiel) then
            write (NFSOUT,7701) iel,NRHS,NDOFM,NFW,KFW,LFW,NE,mke,mkf
 7701        format(1x,'FRWCP: iel,NRHS,NDOFM,NFW,KFW,LFW,NE,mke,mkf',
     .              /,10x,10i8)
         endif
         if (IPFSLH .eq. 1 .or. IPFSLH.eq.negiel) then
            write (NFSOUT,7702) (Elem(i),i=1,mke)
 7702      format(1x,'FRWCP: Elem: Element LHS',/,(10x,1p,10e12.5))
         endif
         if (IPFSRH .eq. 1 .or. IPFSRH.eq.negiel) then
            write (NFSOUT,7703) ((Elem(mke+i),i=1,NDOFM),j=1,NRHS)
 7703      format(1x,'FRWCP: Elem: Element RHS',/,(10x,1p,10e12.5))
         endif
c
c assemble the Element rhs & lhs into the front
c ---------------------------------------------
         if (IASSEM .ne. 0) then
cwb >
c to save having to do abs() repeatedly, do it now
c  **note: sub.dest is recalled in sub.bckwrd so this wont screw it up
c          BUT the negative amdest info is now lost herein
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
            do 133 id = 1,NDOFM
               Amdest(id) = dabs(Amdest(id))
  133       continue
cwb <
c
c symmetric
c
            if (ISYM.eq.1 .or. ISYM.eq.4) then
               call symasm (Amdest, Elem, Frnt)
c             --------------------------------
c unsymmetric
c
            else
               call unsasm (Amdest, Elem, Frnt)
c              --------------------------------
            endif
cwb >
            if (IERR .ne. 0) return
cwb <
c
            if(NRHS .ne. 0) call semrhs (Amdest, Elem(mke+1),
     .                                           Frnt(MKF+1))
c                           ---------------------------------
c debug print
c
            if (IPFSLF .eq. 1 .or. IPFSLF.eq.negiel) then
               lenf = (NFW*(NFW+1))/2
               if (ISYM.eq.2 .or. ISYM.eq.3) lenf = NFW*NFW
               write (NFSOUT,7704) (Frnt(i),i=1,lenf)
 7704           format(1x,'FRWCP: Frnt: Assembled LHS',/,
     .                 (10x,1p,10e12.5))
            endif
            if (IPFSRF .eq. 1 .or. IPFSRF.eq.negiel) then
               write (NFSOUT,7705) ((Frnt(MKF+i),i=1,NFW),j=1,NRHS)
 7705           format(1x,'FRWCP: Frnt: Assembled RHS',/,
     .                 (10x,1p,10e12.5))
            endif
c
         endif
c
c eliminate those equations possible from the current front
c **********************************************************
c
         do 150 ie = 1,ne
c        ****************
c
c check if the equation Buffers are going to get full
c   and dump if necessary
c -------------------------------------------------
c *note: Ubuf runs forward into the workspace Buff(1-->MBUF)
c        Lbuf runs backward into the workspace Buff(1<--MBUF)
c
c
            newu = IU + NFW + NRHS - 1
cwb >
cwb             newl = IL + 1 - NFW * ISYM/3
            newl = IL + 1 - NFW * iunsr
cwb <
c
            if (newu .ge. newl) then
c
c increment the record counter
c
               IFU = IFU + 1
c
c dump the Buffer
c
               call zdirio('U', 'WRITE', IFU, IU-1, Buf, jerr)
c              -----------------------------------------------
               IDUMP = 1
c
c error out
c
               IERR = 10*jerr
               if(jerr .ne. 0) return
c
c set the pointer to the beginning of the equation Buffer
c
               iu = 1
c
c for unsymmetric w/ resolution
c
               if (ISYM .eq. 3) then
                  IFL = IFL + 1
                  call zdirio ('L','WRITE',IFL,MBUF-IL,Buf(IL+1),jerr)
c                 ---------------------------------------------------
                  IL = MBUF
                  IERR = 10*jerr
                  if (jerr .ne. 0) return
               endif
c
c debug print
c
               if (IPFSZD .eq. 1) then
                  write (NFSOUT,7706)  IFU,IFL,IL,IU,MBUF
 7706              format(1x,'FRWCP: IFU,IFL,IL,IU,MBUF',5i8)
               endif
c
c debug print
c
               if (IPFSBF .eq. 1) then
                  write (NFSOUT,7711) (Buf(i),i=1,40)
 7711              format(1x,'FRWCP: First 40 items in Buf:',/,
     .                     (10x,1p,10e12.5))
c
                  if (ISYM .eq. 3) write (NFSOUT,7712)
     .                             (Buf(i),i=IL+1,il+40)
 7712              format(1x,'FRWCP: First 40 items in Buf(IL+1):',/,
     .                    (10x,1p,10e12.5))
c
               endif
c
            endif
c
c eliminate this equation lhs from the current front
c  ----------------------------------------------
            m = IU
c
c symmetric
c
            if (ISYM.eq.1 .or. ISYM.eq.4) then
               call symelm (iel, Andest(ie), Frnt, Buf(IU))
c              --------------------------------------------
c unsymmetric
c
            else
               call unselm (iel, Andest(ie), Frnt, Buf(IU))
c              --------------------------------------------
            endif
c
c debug print
c
            if (IPFSLE .eq. 1 .or. IPFSLE.eq.negiel) then
               lenf = (NFW*(NFW+1))/2
               if (ISYM.eq.2 .or. ISYM.eq.3) lenf = NFW*NFW
               write (NFSOUT,7707) iel,Andest(ie),(Frnt(i),i=1,lenf)
 7707           format(1x,'FRWCP: iel,Andest(ie):',i8,1p,e12.4,
     .                  /,10x,'Frnt: Eliminated LHS',/,
     .                  (10x,1p,10e12.5))
            endif
c
c increment the equation Buffer pointer
c ------------------------------------
c
            IU = IU + NRHS + NFW
c
            if (IERR .eq. 2) then
               write (NFSOUT,7000) iel
 7000           format(2(/), 5x,'ERROR IN FRWCP:',
     .                ' ZERO PIVOT IN ELEMENT:',i4)
               return
            endif
c
c eliminate this equation rhs from the current front
c  ----------------------------------------------
            if(NRHS  .ne.  0) then
c
c symmetric
c
               if (ISYM.eq.1 .or. ISYM.eq.4) then
                  call elmrhs (Andest(ie),   1, Frnt(MKF+1), Buf(m),
     .                         Buf(m+NFW))
c                 ---------------------------------------------------
c unsymmetric
c
               else
                  call elmrhs (Andest(ie), KFW, Frnt(MKF+1), Frnt(NFW),
     .                         Buf(m+NFW))
c                 -----------------------------------------------------
               endif
c
c debug print
c
               if (IPFSRE .eq. 1 .or. IPFSRE.eq.negiel) then
                  write (NFSOUT,7708) iel,Andest(ie),
     .                                ((Frnt(MKF+i),i=1,NFW),j=1,NRHS)
 7708              format(1x,'FRWCP:  iel,Andest(ie):',i8,1p,e12.4,
     .                     /,10x,'Frnt: Eliminated RHS',/,
     .                     (10x,1p,10e12.5))
               endif
c
            endif
c
c for unsymmetric w/ resolution
c   we must copy the the k(i,m)/k(m,m) terms into the Lbuf Buffer space
c   it may be dumped to diskfile for resolution
c   (see sub.unselm and sub.elmrhs for further explanation)
c   *note: Ubuf runs forward into the workspace Buff(1-->MBUF)
c          Lbuf runs backward into the workspace Buff(1<--MBUF)
c
            if (ISYM .eq. 3) then
               m = NFW
               n = NFW-1
               n = max0(n,1)
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c
               do 110 j = 1,n
                  Buf(IL) = Frnt(m)
                  IL = IL - 1
                  m = m + KFW
  110          continue
c
            endif
c
c reduce the front by the equation we just reduced
c
            NFW = NFW - 1
c           -------------
  150    continue
c ***************
c store the final length of the current front (LFW) into Alelm
c
         LFW = NFW
         Alelm(iel) = LFW
c
c for unsymmetric w/ equations being eliminated
c   we must copy the the k(i,m)/k(m,m) terms forward
c   (see sub.unselm and sub.elmrhs for further explanation)
c
c
         if ((ISYM.eq.2 .or. ISYM.eq.3) .and. ne.ne.0) then
c
c n = the beginning frontwith
c m = the final frontwidth
c
            n = KFW
            m = NFW + 1
c
            do 170 i = 2,NFW
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c
               do 160 j = 1,NFW
                  Frnt(m) = Frnt(n+j)
                  m = m + 1
  160          continue
c
               n = n + KFW
  170       continue
c
         endif
c
  200 continue
c ************
c
c dump the final Buffers to diskfiles
c  (or the only Buffers, full or not, for resolution)
c---------------------------------------------------------------
c
      IFU = IFU + 1
cwb >
cwb       call zdirio ( 'u', 'write', IFU, iu-1, Buf, jerr)
cwb c     --------------------------------------------------
cwb
      jerr = 0
c
c if we are not going to do a resolution (ISYM=3 OR 4)
c   and we have not had to dump any Buffers,  then skip this write
c
      if (ISYM.eq.4 .or. IFU.gt.1 .or. ISYM.eq.3) then
         IDUMP = 1
         call zdirio ('U', 'WRITE', IFU, IU-1, Buf, jerr)
c       -------------------------------------------------
c
      else
         LENU = IU - 1
      endif
cwb <
      IERR = 10*jerr
      if(jerr .ne. 0) return
c
      if (ISYM .eq. 3) then
         IFL = IFL+1
         call zdirio ('L', 'WRITE', IFL, MBUF-IL, Buf(IL+1), jerr)
c        ----------------------------------------------------------
         IERR=10*jerr
         if(jerr .ne. 0) return
      endif
c
c debug print
c
      if (IPFSBF .eq. 1) then
         write (NFSOUT,7709) (Buf(i),i=1,40)
 7709     format(1x,'FRWCP: First 40 items in Buf: ',/,
     .                     (10x,1p,10e12.5))
c
         if (ISYM .eq. 3) write (NFSOUT,7710) (Buf(i),i=il+1,il+40)
 7710                   format(1x,'FRWCP: First 40 items in Buf(il+1):',
     .                         /,(10x,1p,10e12.5))
c
      endif
c
c
      end subroutine frwcp
