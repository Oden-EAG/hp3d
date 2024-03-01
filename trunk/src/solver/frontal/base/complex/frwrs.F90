!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:  Forward elimination of rhs's only
!  calls solin1 for dest. vectors
!  calls solin2 for Element rhs's
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
!
! W  Alelm   length of the current front as we process each Elem
! W  Aldest  nodal destination vectors
! W  Amdest  dof destination vectors
! W  Andest  destination in the front where the remaining values
! W          in the front will transfer to as each exiting dof is elimin
! W  Elem    Element lhs + rhs
! W  Frnt    equations in the front
! W  Bbuf    Buffer to contain the new eliminated rhs
! W  Buf     eliminated equations
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
      subroutine frwrs (Alelm, Aldest, Amdest, Andest, &
                        Elem, Frnt, Bbuf, Buf)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      complex(8) :: Aldest(*),Amdest(*),Andest(*),Elem(*),Frnt(*), &
                    Bbuf(*),Buf(*),Alelm(*)
!
      integer :: i,id,ie,iel,ifg,igot,ill,inc
      integer :: irout4,irout5,is,j,jerr,md,n,ne,negiel,numdes
!
      IB = 1
!cwb >
! In order to allow resolution non-sequentially we need to allow IFU to
!  counted herein, not picked up off the common block
!
!cwb       ifui = IFU
!cwb       ifli = IFL
!cwb       IFU = 0
      if (isym .eq. 4) IFU = 0
!cwb <
!
      IFL = 0
      IFB = 0
!
      NFW = 0
      LFW = 0
!cwb >
      IDUMP = 1
      LENU = 0
      NRHSF = NRHS
!cwb <
!
! for symmetric case: we move forward thru the Ubuf buffer to pick up [k
! for unsymmetric case: we move backward thru the Lbuf buffer to pick up
!
      if (ISYM .eq. 4) then
         is = MBUF
         ill = 0
         inc = 1
      else
         is = 1
         ill = MBUF
         inc = -1
      endif




!
! loop over the Elements
! **********************
      do 200 iel  = 1,NUMELM
!     ***********************




!
         negiel = -iel
!
! pull up the Elements nodal destination vectors (Aldest)
!
         call solin1 (iel, numdes, Aldest)
!        ---------------------------------
         if (numdes.eq.0) go to 200
!
! convert the nodal destination vectors to dof destination vectors (Amde
!  also determine the number of dof to eliminate from the front (ne)
!  and the transfer loc for active dof in the front after eliminating (A
!
         call dest (numdes, Aldest, ne, Amdest, Andest)
!        -----------------------------------------------
! check that the frontwidth has not reduced, if so, dont let it
!  this is necessary especially for indexing for unsymmetric k
!
         if (LFW .gt. NFW) NFW  =  LFW
!
! zero out for any new dof in the front
!
         call zeror (Frnt)
!        ----------------
!
! build the Element rhs
!----------------------------
! pointers into Elem(): (1: points to rhs)
! flag for solin2 (ifg=1: get rhs ; =2: get lhs as well)
!
         ifg = 1
         md = int(Amdest(1))
         call solin2 (iel, ifg, NDOFM, NRHS, md, Elem, Elem)
!        -------------------------------------------------------
! debug print
!
         if (IPFSST .eq. 1 .or. IPFSST.eq.negiel) then
            write (NFSOUT,7701) iel,NRHS,NDOFM,NFW,LFW,NE
 7701        format(1x,'FRWRS: iel,NRHS,NDOFM,NFW,LFW,NE', &
                    /,10x,10i8)
         endif
         if (IPFSRH .eq. 1 .or. IPFSRH.eq.negiel) then
            write (NFSOUT,7703) ((Elem(i),i=1,NDOFM),j=1,NRHS)
 7703      format (1x,'FRWRS: Elem: Element RHS',/,(10x,1p,10e12.5))
         endif
!
! assemble the Element rhs & lhs into the front
! ---------------------------------------------
         if (IASSEM .ne. 0) then
!cwb >
! to save having to do abs() repeatedly, do it now
!  **note: sub.dest is recalled in sub.bckwrd so this wont screw it up
!          BUT the negative amdest info is now lost herein
!
!
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!


            do 133 id = 1,NDOFM
cwr10.07.99
ccccc          Amdest(id) = Cdabs(Amdest(id))
               Amdest(id) =   abs(Amdest(id))
  133       continue



!cwb <
!
            call semrhs (Amdest, Elem, Frnt)
!           -----------------------------
! debug print
!
            if (IPFSRF .eq. 1 .or. IPFSRF.eq.negiel) then
               write (NFSOUT,7702) ((Frnt(i),i=1,NFW),j=1,NRHS)
 7702         format(1x,'FRWRS: Frnt: Assembled RHS',/,(10x,1p,10e12.5))
            endif
!
         endif
!
! eliminate those equations possible from the current front
! **********************************************************
!
         do 150 ie = 1,ne
!        ****************

!
! check for rhs Buffer being full and dump if necessary
! -------------------------------------------------------
!
            n  =  IB + NRHS - 1



!
            if (n .gt. MB) then
!
! increment the record counter
!
               IFB = IFB+1
!
! dump the Buffer

!
            if(IDUMPWR.eq.1)then
              write(nfsout,*) &
              'FRWRS: TOO LITTLE SPACE TO AVOID DUMPS'
              stop 'frwrs 1'
            else
               call zdirio ('B', 'WRITE' ,IFB, IB-1, Bbuf, jerr)
!              ----------------------------------------------
            endif
! error out
!
               IERR = 10*jerr
               if (jerr.ne.0) return
!
! set the pointer to the beginning of the rhs Buffer
!
               IB  =  1
!
! debug print
!
               if (IPFSBF .eq. 1) then
                  write (NFSOUT,7704) (Bbuf(i),i=1,40)
 7704              format(1x,'FRWCP: First 40 items in Bbuf:',/, &
                           (10x,1p,10e12.5))
               endif
!
            endif
!
! read back the eliminated lhs equations from Buffers; for symmetric cas
! ----------------------------------------------------------------------
            igot = 0
            if (ISYM .eq. 4) then


               if (is .gt. ill) then
                  IFU = IFU + 1

                  if(IDUMPWR.eq.1)then
                    if(irout4.ne.100)then
                      irout4=100
                      write(nfsout,*)'FRWRS: PRETEND READ U'
                    endif
                    call zdiriodum('U', 'READ', IFU, ill, Buf, jerr)
                  else
                    call zdirio   ('U', 'READ', IFU, ill, Buf, jerr)
!                   ----------------------------------------------
                  endif


                  IERR = 10*jerr
                  is = 1
                  if(jerr.ne.0) return
                  igot = 1
               endif

!
            else
!
! read back the eliminated lhs equations from Buffers; for unsymmetric c
!-----------------------------------------------------------------------
! **Note: for unsymmetric we use both Ubuf and Lbuf for holding the
!          eliminated lhs equations, where Lbuf holds the form of the el
!          used for eliminating the rhs, and Ubuf holds the form of the
!          used for back substitution
!
               if (is .lt. ill) then
                  IFL = IFL + 1


                  call zdirio ('L', 'READ', IFL, ill, Buf, jerr)
!                  ---------------------------------------------
                  IERR = 10*jerr
                  is = ill
                  ill = 1
                  if (jerr.ne.0) return
                  igot = 1
               endif
            endif
!
! debug print
!
           if (IPFSBF .eq. 1 .and. igot.eq.1) then
              write (NFSOUT,7705) (Buf(i),i=1,40)
 7705          format(1x,'FRWCP: First 40 items in Buf: ',/, &
                           (10x,1p,10e12.5))
           endif
!
!
! eliminate this equation rhs from the current front
!  ----------------------------------------------
!
            call elmrhs (Andest(ie), inc, Frnt, Buf(is), Bbuf(ib))
!           -------------------------------------------------------
! debug print
!
             if (IPFSRE .eq. 1 .or. IPFSRE.eq.negiel) then
                write (NFSOUT,7706) iel,Andest(ie), &
                                      ((Frnt(i),i=1,NFW),j=1,NRHS)
 7706            format(1x,'FRWRS:  iel,Andest(ie):',i8,1p,e12.4, &
                         /,10x,'Frnt: Eliminated RHS',/, &
                         (10x,1p,10e12.5))
             endif
!
! increment the rhs Buffer pointer
!
            IB  =  IB + NRHS
!
! increment the lhs equation Buffer pointer
!
            if (ISYM .eq. 4) is = is + NFW + NRHSF
            if (ISYM .eq. 3) is = is - NFW + 1
!
! reduce the front by the equation we just reduced
!
            NFW  =  NFW - 1
!
  150    continue
! ***************
! store the final length of the current front (LFW) into Alelm
!
         LFW  =  NFW
         Alelm(iel)  =  LFW
!
  200 continue
! *************
!
! dump the final rhs Buffers to diskfiles
!---------------------------------------------------------------
!
      IFB = IFB+1
!
      if(IDUMPWR.eq.1)then
        if(irout5.ne.100)then
          write(nfsout,*)'FRWRS: PRETEND DUMPING B'
          irout5=100
        endif
        call zdiriodum('B', 'WRITE', IFB, IB-1, Bbuf, jerr)
      else
        call zdirio  ('B', 'WRITE', IFB, IB-1, Bbuf, jerr)
!       ---------------------------------------------
      endif
!
      IERR = 10*jerr
!
      IB = 1
!cwb >
!cwb       IFU = ifui
!cwb       IFL = ifli
!cwb
!
! debug print
!
      if (IPFSBF .eq. 1) then
         write (NFSOUT,7707) (Bbuf(i),i=1,40)
 7707     format(1x,'FRWRS: First 40 items in Bbuf: ',/, &
                           (10x,1p,10e12.5))
      endif
!
!
      end subroutine frwrs
