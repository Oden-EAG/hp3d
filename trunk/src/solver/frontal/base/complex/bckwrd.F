      subroutine bckwrd (Alelm, Aldest, Amdest,
     .                   Andest, Elem, Frnt, Bbuf, Ubuf)
c
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c     Backsubstitution
c     -----------------
c     calls SOLIN1 for dest. vectors
c     passes Elemental solutions to SOLOUT
**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c W  Alelm   final length of the current front as we process each Elem
c W  Aldest  nodal destination vectors
c W  Amdest  dof destination vectors
c W  Andest  Destination in the front where the remaining values
c            in the front will transfer to as each exiting dof is elimin
c
c W  Elem    Temp storage for the holding the Element solution
c             (ie: used to transfer x from the solution front to the
c                  local Element numbering based on the dest vecs)
c
c W  Frnt    Temp storage for the solutions
c            Note: on entry, this contains the previous values in the fr
c                   the values at the end of sub.frwcp or sub.frwrs
c                   but this info isnt used for anything
c
c W  Ubuf    The eliminated lhs buffered equations
c W  Bbuf    The eliminated rhs buffered equations
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
c
      use surfsc1
      use surfsc2
#include "syscomz.blk"
      dimension Aldest(*), Amdest(*), Andest(*), Elem(*), Frnt(*),
     .          Bbuf(*), Ubuf(*), Alelm(*)


c  ...test......test......test......test......test......test......test...
ccc      integer :: ndest1(100)
c  ...test......test......test......test......test......test......test...

c
#if DEBUG_MODE
      integer :: iprint=0
#endif
c
      iuu = 1
      jel = NUMELM + 1
c
      IB = 0
      ifui = IFU
      ifbi = IFB
c
      IFU = IFU + 1
      IFB = IFB + 1

c
c loop over the Elements (we are actually looping thru them backwards)
c **********************
      do 200 iel = 1,NUMELM
c     *******************
c set to loop backwards thru the Elems

c
         jel = jel - 1
c
         negiel = -jel
c
c pull up the Elements nodal destination vectors (Aldest)
c

         call solin1 (jel, numdes, Aldest)
c        ----------------------------------
#if DEBUG_MODE
         if (iprint.eq.1) then
           write(NFSOUT,*) 'IN BCKWRD BEFORE CALL TO dest= '
           write(NFSOUT,*) 'jel,numdes = ',jel,numdes
           write(NFSOUT,*) 'Aldest = ',(Aldest(ii),ii=1,numdes)
         endif
#endif

c  ...test......test......test......test......test......test......test...
c  ...test......test......test......test......test......test......test...
c         read(31,*)iel1
c         read(31,*)numdes1
c         read(31,*)(ndest1(i),i=1,numdes1)
c         if(numdes.ne.numdes1)then
c           write(*,*)'BCKWRD: iel,numdes,numdes1=',iel,numdes,numdes1
c           stop
c         endif
c         do i=1,numdes
c            if(ndest1(i) .ne. nint(dreal(aldest(i))))then
c              write(*,*)'BCKWRD:IEL,I=',iel,i
c              write(*,*)ndest1(i),nint(dreal(aldest(i)))
c              stop
c            endif
c            enddo
c  ...test......test......test......test......test......test......test...
c  ...test......test......test......test......test......test......test...


c
         if (numdes .eq. 0) go to 200
c
c convert the nodal destination vectors to dof destination vectors (Amde
c  also determine the number of dof to eliminate from the front (ne)
c  and the transfer loc for active dof in the front after eliminating (A
c
         call dest (numdes, Aldest, ne ,Amdest, Andest)
c        -----------------------------------------------
#if DEBUG_MODE
         if (iprint.eq.1) then
           write(NFSOUT,*) 'IN BCKWRD NDOFM AFTER CALL TO dest= ',NDOFM
         endif
#endif
c
c pick up final frontwidth (after elim) for the next Element
c
         if (jel .ne. 1) then
            LFW = Alelm(jel-1)
            if (LFW .gt. NFW) NFW = LFW
         endif
c
         NFW = NFW - ne + 1
         j = ne + 1
c
         if (IPFSST .eq. 1 .or. IPFSST.eq.negiel) then
            write (NFSOUT,7701) jel,NRHS,NDOFM,NFW,LFW,NE
 7701        format(1x,'BCKWRD: jel,NRHS,NDOFM,NFW,LFW,NE',
     .              /,10x,10i8)
         endif


c
c backsubstitute those equations possible from the current front
c **********************************************************
         do 150 ie = 1,ne
c        ****************
            j = j - 1



c
c***note:  we loop backwards thru the buffered equations (iuu runs backw
c          also nfw grows as we proceed  (1,2,3,...)
c          we also run backwards thru Andest()
c
c if we need to read up a buffer,
c -------------------------------
c   (and in fact we had to write one previously)
c
            if (iuu .le. 1 ) then
c
c if we act ually dumped buffers, then read them up
c  Note: for resolution w/ unsymmetric we dont know IFU apriori
c         and zdirio passes it back
c
cwb >
               if (IDUMP.eq.1 .or. IRESOL.eq.1) then
cwb <
                 IFU = IFU - 1

                 if(IDUMPWR.eq.1)then
                   if(irout4.ne.100)then
                     write(NFSOUT,*)'BCKWRD: PRETEND READ U'
                     irout4=100
                   endif
                     call zdiriodum ('U', 'READ', IFU, ill, Ubuf, jerr)
                 else
                     call zdirio    ('U', 'READ', IFU, ill, Ubuf, jerr)
                 endif

cwb >
c for unsymmetric during resolution we didnt know IFU apriori
c  and it is sent back from zdirio
c
                 if (ISYM.eq.3 .and. IRESOL.eq.1)  ifui = IFU
cwb <
                 IERR = 10*jerr
                 if (jerr .ne. 0) return
c
                 iuu = ill + 1
cwb >
c if we didnt dump buffers, just set the read flags
c
              else
                 IFU = IFU - 1
                 iuu = LENU + 1
              endif
cwb <
            endif
c
c calculate the solution for this dof
c -----------------------------------
c iuu - pointer into the lhs eq buffer
c u  = { ..... |k1,k2,k3, | k1,k2, | k1 }
c                                  <=== iuu:  runs backwards up the list
c
            iuuold = iuu
            iuu = iuu - NFW - NRHSF
c
c debug print
c
            if (IPFSBK .eq. 1 .or. IPFSBK.eq.negiel) then
               write (NFSOUT,7702) jel,ie,Andest(j)
 7702           format(1x,'BCKWRD: jel,ie,Andest(j)',2i8,1p,e12.5)
               write(NFSOUT,7703) (Ubuf(i),i=iuu,iuuold)
 7703           format(1x,'BCKWRD: (Ubuf(i),i=iuu,iuuold)',/,
     .                (10x,1p,10e12.5))
            endif
c
            if (IRESOL .ne. 1) then
c
c calculate the solution (no resolution)
c
               n = iuu + NFW

               call elmsol (Andest(j), Ubuf(iuu), Ubuf(n), Frnt)
c              -------------------------------------------------
c debug print
c
               if (IPFSBK .eq. 1 .or. IPFSBK.eq.negiel) then
                  write(NFSOUT,7704) (Ubuf(n+i),i=1,NRHS)
 7704              format(1x,'BCKWRD: (Ubuf(n+i),i=1,NRHS)',
     .                    (/,10x,1p,10e12.5))
               endif
c
            else
c
c if this is a resolution, we must read back in the recalculated
c  elimination buffers for the rhs (from sub.frwrs)
c
               if (IB .le. 0) then
                  IFB = IFB - 1

              if(IDUMPWR.eq.1)then
                if(irout5.ne.100)then
                  write(NFSOUT,*)'BCKWRD: PRETEND READ B'
                  irout5=100
                endif
                  call zdiriodum('B', 'READ', IFB, ill, Bbuf, jerr)
              else
                  call zdirio   ('B', 'READ', IFB, ill, Bbuf, jerr)
              endif

c error out
c
                  IERR = 10*jerr
                  if (jerr .ne .0) return
c
                  IB = ill - NRHS + 1
               endif
c
c calculate the solution (resolution)
c
               call elmsol (Andest(j), Ubuf(iuu), Bbuf(IB), Frnt)
c              ---------------------------------------------
c debug print
c
               if (IPFSBK .eq. 1 .or. IPFSBK.eq.negiel) then
                  write(NFSOUT,7705) (Bbuf(IB+i),i=1,NRHS)
 7705             format(1x,'BCKWRD: (Bbuf(IB+i),i=1,NRHS)',/,
     .                    (10x,1p,10e12.5))
               endif
c
               IB = IB - NRHS
c
            endif
c
            NFW = NFW + 1
  150    continue
c
c pass the solution to the driving program for storage
c ----------------------------------------------------
c
         if (IASSEM .eq. 0) then
c
            call solout (jel, NDOFM, NRHS, Amdest, Frnt)
c           ---------------------------------------------
         else
c
cwb >
c transfer x from the solution front to the
c  local Element numbering based on the dest vecs
cwb* **note: 1) we may want to inline locr below for greater speed
c            2) we also should switch the loop parameters below
c
            k = 1
            do 180 j = 1, NRHS
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
               do 160 i = 1,NDOFM
                  md = Amdest(i)
                  md = abs(md)
                  l = (j-1)*MFW + md
                  Elem(k) = Frnt(l)
                  k = k + 1
  160          continue
c
  180       continue
c
cwb             do 180 i = 1,NDOFM
cwb                k = 0
cwb c
cwb cvd$ select (vector)
cwb c
cwb                do 160 j = 1,NRHS
cwb                   md = Amdest(i)
cwb                   l = locr(j,md)
cwb                   Elem(k+i) = Frnt(l)
cwb                   k = k + NDOFM
cwb   160          continue
cwb c
cwb   180       continue
cwb <
c
#if DEBUG_MODE
            if (iprint.eq.1) then
              write(NFSOUT,*) 'IN BCKWRD jel,NDOFM = ',jel,NDOFM
            endif
#endif
            call solout (jel, NDOFM, NRHS, Amdest, Elem)
c           ---------------------------------------------
c debug print
c
               if (IPFSXX .eq. 1 .or. IPFSXX.eq.negiel) then
                  lenx = NDOFM*NRHS
                  write(NFSOUT,7706) (Elem(i),i=1,lenx)
 7706              format(1x,'BCKWRD: (Elem(i),i=1,lenx)',/,
     .                    (10x,1p,10e12.5))
               endif
c
         endif
c
  200 continue
c ************
c
      IFU = ifui
      IFB = ifbi
c
c
      end subroutine bckwrd
