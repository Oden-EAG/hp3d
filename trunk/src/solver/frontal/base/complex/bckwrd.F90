!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
!     Backsubstitution
!     -----------------
!     calls SOLIN1 for dest. vectors
!     passes Elemental solutions to SOLOUT
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! W  Alelm   final length of the current front as we process each Elem
! W  Aldest  nodal destination vectors
! W  Amdest  dof destination vectors
! W  Andest  Destination in the front where the remaining values
!            in the front will transfer to as each exiting dof is elimin
!
! W  Elem    Temp storage for the holding the Element solution
!             (ie: used to transfer x from the solution front to the
!                  local Element numbering based on the dest vecs)
!
! W  Frnt    Temp storage for the solutions
!            Note: on entry, this contains the previous values in the fr
!                   the values at the end of sub.frwcp or sub.frwrs
!                   but this info isnt used for anything
!
! W  Ubuf    The eliminated lhs buffered equations
! W  Bbuf    The eliminated rhs buffered equations
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
      subroutine bckwrd (Alelm, Aldest, Amdest, &
                         Andest, Elem, Frnt, Bbuf, Ubuf)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      complex(8) :: Aldest(*), Amdest(*), Andest(*), Elem(*), &
                    Frnt(*), Bbuf(*), Ubuf(*), Alelm(*)
!
      integer :: iuu,jel,i,ie,iel,ifbi,ifui,ill
      integer :: irout4,irout5,iuuold,j,jerr,k,l,lenx
      integer :: md,n,ne,negiel,numdes
!
!  ...test......test......test......test......test......test......test...
ccc      integer :: ndest1(100)
!  ...test......test......test......test......test......test......test...

!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
      iuu = 1
      jel = NUMELM + 1
!
      IB = 0
      ifui = IFU
      ifbi = IFB
!
      IFU = IFU + 1
      IFB = IFB + 1

!
! loop over the Elements (we are actually looping thru them backwards)
! **********************
      do 200 iel = 1,NUMELM
!     *******************
! set to loop backwards thru the Elems

!
         jel = jel - 1
!
         negiel = -jel
!
! pull up the Elements nodal destination vectors (Aldest)
!

         call solin1 (jel, numdes, Aldest)
!        ----------------------------------
#if HP3D_DEBUG
         if (iprint.eq.1) then
           write(NFSOUT,*) 'IN BCKWRD BEFORE CALL TO dest= '
           write(NFSOUT,*) 'jel,numdes = ',jel,numdes
           write(NFSOUT,*) 'Aldest = ',(Aldest(ii),ii=1,numdes)
         endif
#endif

!  ...test......test......test......test......test......test......test...
!  ...test......test......test......test......test......test......test...
!         read(31,*)iel1
!         read(31,*)numdes1
!         read(31,*)(ndest1(i),i=1,numdes1)
!         if(numdes.ne.numdes1)then
!           write(*,*)'BCKWRD: iel,numdes,numdes1=',iel,numdes,numdes1
!           stop
!         endif
!         do i=1,numdes
!            if(ndest1(i) .ne. nint(dreal(aldest(i))))then
!              write(*,*)'BCKWRD:IEL,I=',iel,i
!              write(*,*)ndest1(i),nint(dreal(aldest(i)))
!              stop
!            endif
!            enddo
!  ...test......test......test......test......test......test......test...
!  ...test......test......test......test......test......test......test...


!
         if (numdes .eq. 0) go to 200
!
! convert the nodal destination vectors to dof destination vectors (Amde
!  also determine the number of dof to eliminate from the front (ne)
!  and the transfer loc for active dof in the front after eliminating (A
!
         call dest (numdes, Aldest, ne ,Amdest, Andest)
!        -----------------------------------------------
#if HP3D_DEBUG
         if (iprint.eq.1) then
           write(NFSOUT,*) 'IN BCKWRD NDOFM AFTER CALL TO dest= ',NDOFM
         endif
#endif
!
! pick up final frontwidth (after elim) for the next Element
!
         if (jel .ne. 1) then
            LFW = int(Alelm(jel-1))
            if (LFW .gt. NFW) NFW = LFW
         endif
!
         NFW = NFW - ne + 1
         j = ne + 1
!
         if (IPFSST .eq. 1 .or. IPFSST.eq.negiel) then
            write (NFSOUT,7701) jel,NRHS,NDOFM,NFW,LFW,NE
 7701        format(1x,'BCKWRD: jel,NRHS,NDOFM,NFW,LFW,NE', &
                    /,10x,10i8)
         endif


!
! backsubstitute those equations possible from the current front
! **********************************************************
         do 150 ie = 1,ne
!        ****************
            j = j - 1



!
!***note:  we loop backwards thru the buffered equations (iuu runs backw
!          also nfw grows as we proceed  (1,2,3,...)
!          we also run backwards thru Andest()
!
! if we need to read up a buffer,
! -------------------------------
!   (and in fact we had to write one previously)
!
            if (iuu .le. 1 ) then
!
! if we act ually dumped buffers, then read them up
!  Note: for resolution w/ unsymmetric we dont know IFU apriori
!         and zdirio passes it back
!
!cwb >
               if (IDUMP.eq.1 .or. IRESOL.eq.1) then
!cwb <
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

!cwb >
! for unsymmetric during resolution we didnt know IFU apriori
!  and it is sent back from zdirio
!
                 if (ISYM.eq.3 .and. IRESOL.eq.1)  ifui = IFU
!cwb <
                 IERR = 10*jerr
                 if (jerr .ne. 0) return
!
                 iuu = ill + 1
!cwb >
! if we didnt dump buffers, just set the read flags
!
              else
                 IFU = IFU - 1
                 iuu = LENU + 1
              endif
!cwb <
            endif
!
! calculate the solution for this dof
! -----------------------------------
! iuu - pointer into the lhs eq buffer
! u  = { ..... |k1,k2,k3, | k1,k2, | k1 }
!                                  <=== iuu:  runs backwards up the list
!
            iuuold = iuu
            iuu = iuu - NFW - NRHSF
!
! debug print
!
            if (IPFSBK .eq. 1 .or. IPFSBK.eq.negiel) then
               write (NFSOUT,7702) jel,ie,Andest(j)
 7702           format(1x,'BCKWRD: jel,ie,Andest(j)',2i8,1p,e12.5)
               write(NFSOUT,7703) (Ubuf(i),i=iuu,iuuold)
 7703           format(1x,'BCKWRD: (Ubuf(i),i=iuu,iuuold)',/, &
                      (10x,1p,10e12.5))
            endif
!
            if (IRESOL .ne. 1) then
!
! calculate the solution (no resolution)
!
               n = iuu + NFW

               call elmsol (Andest(j), Ubuf(iuu), Ubuf(n), Frnt)
!              -------------------------------------------------
! debug print
!
               if (IPFSBK .eq. 1 .or. IPFSBK.eq.negiel) then
                  write(NFSOUT,7704) (Ubuf(n+i),i=1,NRHS)
 7704              format(1x,'BCKWRD: (Ubuf(n+i),i=1,NRHS)', &
                          (/,10x,1p,10e12.5))
               endif
!
            else
!
! if this is a resolution, we must read back in the recalculated
!  elimination buffers for the rhs (from sub.frwrs)
!
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

! error out
!
                  IERR = 10*jerr
                  if (jerr .ne .0) return
!
                  IB = ill - NRHS + 1
               endif
!
! calculate the solution (resolution)
!
               call elmsol (Andest(j), Ubuf(iuu), Bbuf(IB), Frnt)
!              ---------------------------------------------
! debug print
!
               if (IPFSBK .eq. 1 .or. IPFSBK.eq.negiel) then
                  write(NFSOUT,7705) (Bbuf(IB+i),i=1,NRHS)
 7705             format(1x,'BCKWRD: (Bbuf(IB+i),i=1,NRHS)',/, &
                          (10x,1p,10e12.5))
               endif
!
               IB = IB - NRHS
!
            endif
!
            NFW = NFW + 1
  150    continue
!
! pass the solution to the driving program for storage
! ----------------------------------------------------
!
         if (IASSEM .eq. 0) then
!
            md = int(Amdest(1))
            call solout (jel, NDOFM, NRHS, md, Frnt)
!           ---------------------------------------------
         else
!
!cwb >
! transfer x from the solution front to the
!  local Element numbering based on the dest vecs
!cwb* **note: 1) we may want to inline locr below for greater speed
!            2) we also should switch the loop parameters below
!
            k = 1
            do 180 j = 1, NRHS
!
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!
               do 160 i = 1,NDOFM
                  md = int(Amdest(i))
                  md = abs(md)
                  l = (j-1)*MFW + md
                  Elem(k) = Frnt(l)
                  k = k + 1
  160          continue
!
  180       continue
!
!cwb             do 180 i = 1,NDOFM
!cwb                k = 0
!cwb !
!cwb !cvd$ select (vector)
!cwb !
!cwb                do 160 j = 1,NRHS
!cwb                   md = int(Amdest(i))
!cwb                   l = locr(j,md)
!cwb                   Elem(k+i) = Frnt(l)
!cwb                   k = k + NDOFM
!cwb   160          continue
!cwb !
!cwb   180       continue
!cwb <
!
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(NFSOUT,*) 'IN BCKWRD jel,NDOFM = ',jel,NDOFM
            endif
#endif
            md = int(Amdest(1))
            call solout (jel, NDOFM, NRHS, md, Elem)
!           ---------------------------------------------
! debug print
!
               if (IPFSXX .eq. 1 .or. IPFSXX.eq.negiel) then
                  lenx = NDOFM*NRHS
                  write(NFSOUT,7706) (Elem(i),i=1,lenx)
 7706              format(1x,'BCKWRD: (Elem(i),i=1,lenx)',/, &
                          (10x,1p,10e12.5))
               endif
!
         endif
!
  200 continue
! ************
!
      IFU = ifui
      IFB = ifbi
!
!
      end subroutine bckwrd
