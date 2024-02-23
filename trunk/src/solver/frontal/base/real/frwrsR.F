c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:  Forward elimination of rhs's only
c  calls solin1 for dest. vectors
c  calls solin2 for Element rhs's
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c
c W  Alelm   length of the current front as we process each Elem
c W  Aldest  nodal destination vectors
c W  Amdest  dof destination vectors
c W  Andest  destination in the front where the remaining values
c W          in the front will transfer to as each exiting dof is elimin
c W  Elem    Element lhs + rhs
c W  Frnt    equations in the front
c W  Bbuf    Buffer to contain the new eliminated rhs
c W  Buf     eliminated equations
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
      subroutine frwrs (Alelm, Aldest, Amdest, Andest,
     .                  Elem, Frnt, Bbuf, Buf)
c
      use surfsc1
      use surfsc2
c
      implicit none
c
      real(8) :: Alelm(*),Aldest(*),Amdest(*),Andest(*)
      real(8) :: Elem(*),Frnt(*),Bbuf(*),Buf(*)
c
      integer :: i,id,ie,iel,ifg,igot,ill,inc,is,j,jerr
      integer :: md,n,ne,negiel,numdes
c
      IB = 1
cwb >
c In order to allow resolution non-sequentially we need to allow IFU to
c  counted herein, not picked up off the common block
c
cwb       ifui = IFU
cwb       ifli = IFL
cwb       IFU = 0
      if (isym .eq. 4) IFU = 0
cwb <
c
      IFL = 0
      IFB = 0
c
      NFW = 0
      LFW = 0
cwb >
      IDUMP = 1
      LENU = 0
      NRHSF = NRHS
cwb <
c
c for symmetric case: we move forward thru the Ubuf buffer to pick up [k
c for unsymmetric case: we move backward thru the Lbuf buffer to pick up
c
      if (ISYM .eq. 4) then
         is = MBUF
         ill = 0
         inc = 1
      else
         is = 1
         ill = MBUF
         inc = -1
      endif
c
c loop over the Elements
c **********************
      do 200 iel  = 1,NUMELM
c     ***********************
c
         negiel = -iel
c
c pull up the Elements nodal destination vectors (Aldest)
c
         call solin1 (iel, numdes, Aldest)
c        ---------------------------------
         if (numdes.eq.0) go to 200
c
c convert the nodal destination vectors to dof destination vectors (Amde
c  also determine the number of dof to eliminate from the front (ne)
c  and the transfer loc for active dof in the front after eliminating (A
c
         call dest (numdes, Aldest, ne, Amdest, Andest)
c        -----------------------------------------------
c check that the frontwidth has not reduced, if so, dont let it
c  this is necessary especially for indexing for unsymmetric k
c
         if (LFW .gt. NFW) NFW  =  LFW
c
c zero out for any new dof in the front
c
         call zeror (Frnt)
c        ----------------
c
c build the Element rhs
c----------------------------
c pointers into Elem(): (1: points to rhs)
c flag for solin2 (ifg=1: get rhs ; =2: get lhs as well)
c
         ifg = 1
         md = int(Amdest(1))
         call solin2 (iel, ifg, NDOFM, NRHS, md, Elem, Elem)
c        -------------------------------------------------------
c debug print
c
         if (IPFSST .eq. 1 .or. IPFSST.eq.negiel) then
            write (NFSOUT,7701) iel,NRHS,NDOFM,NFW,LFW,NE
 7701        format(1x,'FRWRS: iel,NRHS,NDOFM,NFW,LFW,NE',
     .              /,10x,10i8)
         endif
         if (IPFSRH .eq. 1 .or. IPFSRH.eq.negiel) then
            write (NFSOUT,7703) ((Elem(i),i=1,NDOFM),j=1,NRHS)
 7703      format (1x,'FRWRS: Elem: Element RHS',/,(10x,1p,10e12.5))
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
            call semrhs (Amdest, Elem, Frnt)
c           -----------------------------
c debug print
c
            if (IPFSRF .eq. 1 .or. IPFSRF.eq.negiel) then
               write (NFSOUT,7702) ((Frnt(i),i=1,NFW),j=1,NRHS)
 7702         format(1x,'FRWRS: Frnt: Assembled RHS',/,(10x,1p,10e12.5))
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
c check for rhs Buffer being full and dump if necessary
c -------------------------------------------------------
c
            n  =  IB + NRHS - 1
c
            if (n .gt. MB) then
c
c increment the record counter
c
               IFB = IFB+1
c
c dump the Buffer
c
               call zdirio ('B', 'WRITE' ,IFB, IB-1, Bbuf, jerr)
c              ----------------------------------------------
c error out
c
               IERR = 10*jerr
               if (jerr.ne.0) return
c
c set the pointer to the beginning of the rhs Buffer
c
               IB  =  1
c
c debug print
c
               if (IPFSBF .eq. 1) then
                  write (NFSOUT,7704) (Bbuf(i),i=1,40)
 7704              format(1x,'FRWCP: First 40 items in Bbuf:',/,
     .                     (10x,1p,10e12.5))
               endif
c
            endif
c
c read back the eliminated lhs equations from Buffers; for symmetric cas
c ----------------------------------------------------------------------
            igot = 0
            if (ISYM .eq. 4) then
c
               if (is .gt. ill) then
                  IFU = IFU + 1
                  call zdirio ('U', 'READ', IFU, ill, Buf, jerr)
c                  ----------------------------------------------
                  IERR = 10*jerr
                  is = 1
                  if(jerr.ne.0) return
                  igot = 1
               endif
c
            else
c
c read back the eliminated lhs equations from Buffers; for unsymmetric c
c-----------------------------------------------------------------------
c **Note: for unsymmetric we use both Ubuf and Lbuf for holding the
c          eliminated lhs equations, where Lbuf holds the form of the el
c          used for eliminating the rhs, and Ubuf holds the form of the
c          used for back substitution
c
               if (is .lt. ill) then
                  IFL = IFL + 1
                  call zdirio ('L', 'READ', IFL, ill, Buf, jerr)
c                  ---------------------------------------------
                  IERR = 10*jerr
                  is = ill
                  ill = 1
                  if (jerr.ne.0) return
                  igot = 1
               endif
            endif
c
c debug print
c
           if (IPFSBF .eq. 1 .and. igot.eq.1) then
              write (NFSOUT,7705) (Buf(i),i=1,40)
 7705          format(1x,'FRWCP: First 40 items in Buf: ',/,
     .                     (10x,1p,10e12.5))
           endif
c
c
c eliminate this equation rhs from the current front
c  ----------------------------------------------
c
            call elmrhs (Andest(ie), inc, Frnt, Buf(is), Bbuf(ib))
c           -------------------------------------------------------
c debug print
c
             if (IPFSRE .eq. 1 .or. IPFSRE.eq.negiel) then
                write (NFSOUT,7706) iel,Andest(ie),
     .                                ((Frnt(i),i=1,NFW),j=1,NRHS)
 7706            format(1x,'FRWRS:  iel,Andest(ie):',i8,1p,e12.4,
     .                   /,10x,'Frnt: Eliminated RHS',/,
     .                   (10x,1p,10e12.5))
             endif
c
c increment the rhs Buffer pointer
c
            IB  =  IB + NRHS
c
c increment the lhs equation Buffer pointer
c
            if (ISYM .eq. 4) is = is + NFW + NRHSF
            if (ISYM .eq. 3) is = is - NFW + 1
c
c reduce the front by the equation we just reduced
c
            NFW  =  NFW - 1
c
  150    continue
c ***************
c store the final length of the current front (LFW) into Alelm
c
         LFW  =  NFW
         Alelm(iel)  =  LFW
c
  200 continue
c *************
c
c dump the final rhs Buffers to diskfiles
c---------------------------------------------------------------
c
      IFB = IFB+1
c
      call zdirio('B', 'WRITE', IFB, IB-1, Bbuf, jerr)
c     ---------------------------------------------
c
      IERR = 10*jerr
c
      IB = 1
cwb >
cwb       IFU = ifui
cwb       IFL = ifli
cwb
c
c debug print
c
      if (IPFSBF .eq. 1) then
         write (NFSOUT,7707) (Bbuf(i),i=1,40)
 7707     format(1x,'FRWRS: First 40 items in Bbuf: ',/,
     .                     (10x,1p,10e12.5))
      endif
c
c
      end subroutine frwrs
