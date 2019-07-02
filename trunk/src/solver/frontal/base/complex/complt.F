      subroutine complt (Awrk)
c
#include "syscomz.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c Initiate forward elimination of lhs and rhs followed by backsubstituti
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c
c W   Awrk    Workspace used by the solver (it is cut and passed below)
c
c*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
c AUTHOR:                             LATEST REVISION:
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
#include "surfscwr.blk"
c
      dimension Awrk(*)



c
      call zecond(t0)
c
      NNEGP = 0
      NPOSP = 0
      IERR = 0
c
      iprint = 0
c
c set up the storage cut lengths
c--------------------------------
c MDOF = maximum number of dof for any elem
c n = amount of space required to hold the first few temp arrays cut fro
c      the workspace
c
      n = NUMELM + MLDEST + 2*MDOF
c
      if (iprint.eq.1) then
         write(NFSOUT,*) 'in complt.f'
         write(NFSOUT,*) 'mdof,numelm,mldest,n',mdof,numelm,mldest,n
      endif
c MELEM = space required to hold the maximum elem lhs + rhs, during asse
c MKF = space required to hold the lhs equations for the max frontwidth
c
c symmetric
      if(ISYM.eq.1 .or. ISYM.eq.4) then
         MELEM = MDOF*(MDOF+1)/2 + MDOF*NRHS
         MKF = MFW*(MFW+1)/2
c
c unsymmetric
      else
         MELEM = MDOF*(MDOF+NRHS)
         MKF = MFW*MFW
      endif
c
c MFWR = space required to hold the lhs & rhs equations for the max fron
c
      MFWR = MKF + MFW*NRHS
c
c no elem assembly space is needed when the data comes from outside
c  (ie: IASSEM=0)
c
      MELEM = IASSEM*MELEM
c
c MW = total for front and individual elem temp storage
c
      MW = MELEM + MFWR
c
c space which is left over is allocated for buffering equations
c
      MBUF = MA - MW - n
      if (MBUF.lt.1000) then
        write(NFSOUT,*) 'COMPLT:TOO LITTLE BUFFER !! MBUF = ',MBUF
        stop 'complt 1'
      endif
cwb >
c save the value of MBUF & MW to check against during resolution
c *** Note: we cannot let the size of MBUF shrink during resolution
c            because then the stuff that was buffered out will not fit b
c
      if (ISYM .eq. 3) then
         MBUFSV(1,1) = MBUF
         MBUFSV(1,2) = MW
      elseif (ISYM .eq. 4) then
         MBUFSV(2,1) = MBUF
         MBUFSV(2,2) = MW
      endif



cwb <
c
c open direct access buffer files
c ===============================
c open the LHS storage buffer diskfile
c
      if(IDUMPWR.eq.1)then
        write(NFSOUT,*)'COMPLT: PRETEND OPEN U'
        call zdiriodum('U', 'OPEN', 0, MBUF, Awrk, IERR)
      else
        call zdirio ('U', 'OPEN', 0, MBUF, Awrk, IERR)
c       --------------------------------------------
      endif

c
C open L : The additional eliminated LHS buffer for unsymmetric w/ resol
c
      if (ISYM .eq .3) call zdirio ('L', 'OPEN', 0, MBUF, Awrk, IERR)
c                     -------------------------------------------
      IERR = 10*IERR
      if (IERR .ne. 0) go to 9100
c
c debug print
c
      if (IPRSTR .eq. 1) then
         if(ISYM.eq.1 .or. ISYM.eq.4) then
           write(NFSOUT,7000)
         else
           write(NFSOUT,7010)
         endif
         if(ISYM.eq.2 .or. ISYM.eq.1) write(NFSOUT,7020)
         write(NFSOUT,7030) n, MELEM, MFWR, MBUF, MA
c
 7000    format(2(/), 5x,'SYMMETRIC FORWARD ELIMINATION',/)
 7010    format(2(/), 5x,'UNSYMMETRIC FORWARD ELIMINATION',/)
 7020    format(      5x,'RESOLUTION INACTIVATED',/)
 7030    format(     4x,'            OVER HEAD:',i7,/,
     .               4x,'              ELEMENT:',i7,/,
     .               4x,'                FRONT:',i7,/,
     .               4x,'               BUFFER:',i7,/,
     .               4x,'        TOTAL STORAGE:',i7)
      endif
c
c check that the buffer workspace is large enuf
c   it must be able to at least hold the max frontwidth and NRHS rightha
c
      if(MBUF .lt. MFW+NRHS) then
         IERR = 6
         minn = MA + MFW + NRHS - MBUF
         write(NFSOUT,7060) minn
 7060    format(2(/), 5x,'ERROR IN SURFS: NOT ENOUGH ROOM IN BUFFER',
     .             /,12x,'MINIMUM SIZE OF WORK ARRAY  = ',i7)
         go to 9100
      endif
c
c set pointers into the workspace (note: relative pointers)
c ===============================
c set pointer for temp storage of the length of the current front
c   for each elem                                                  (alel
c
      ia1 = 1
c
c set pointer to temp storage for nodal destination vectors       (aldes
c
      ial  =  ia1 + NUMELM
c
c set pointer to temp storage for dof destination vectors         (amdes
c
      iam  =  ial + MLDEST
c
c set pointer to temp storage for the transfer loc for active dof
c   in the front after eliminating                                (andes
c
      ian  =  iam + MDOF
c
c set pointer to temp storage for elem lhs + rhs, during assembly   (ele
c
      iae  =  ian + MDOF
c
c set pointer to temp storage for lhs & rhs equations for front     (frn
c
      iaf  =  iae + MELEM
c
c set pointer to temp storage for equation buffering,
c   this is all remaining workspace                                  (bu
c
      iab  =  iaf + MFWR
c
c debug print
c
      if (IPFSST .eq. 1) then
         write(NFSOUT,7701) NUMELM,MLDEST,MDOF,NRHS,MELEM,MKF,MFWR
 7701     format(1x,'COMPLT: NUMELM,MLDEST,MDOF,NRHS,MELEM,MKF,MFWR',
     .           /,10x,8i10)
         write(NFSOUT,7702) MBUF,MA,MW,ia1,ial,iam,ian,iae,iaf,iab
 7702     format(1x,'COMPLT: MBUF,MA,MW,ia1,ial,iam,ian,iae,iaf,iab',
     .           /,10x,10i10) 
      endif



c
c forward elimination of all elements  (lhs and rhs)
c ====================================
c
      call frwcp (Awrk(ia1), Awrk(ial), Awrk(iam), Awrk(ian),
     .            Awrk(iae), Awrk(iaf), Awrk(iab))
c     --------------------------------------------------------
c
      call zecond(tf)
      dt  =  tf-t0
c
      if(IPRSTR .eq. 1)then
         write(NFSOUT,7080) dt
         iful = IFU*IDUMP + IFL
         if(iful .ne. 0) write(NFSOUT,7100) iful
 7080     format(10x,'TIME IN FORWARD ELIMINATION:',f9.3)
 7100     format(10x,'BUFFER DUMPS:',i4)
      endif
c
      if (IERR .ne. 0) go to 9100
c
      if (NRHS .ne. 0) then


c
c
c backward substitution of all elements
c ====================================
c
         call bckwrd (Awrk(ia1), Awrk(ial), Awrk(iam), Awrk(ian),
     .                Awrk(iae), Awrk(iaf), Awrk(iab), Awrk(iab))
c        --------------------------------------------------------
c
         call zecond(tb)
         dt  =  tb - tf
         if(IPRSTR .eq. 1) write (NFSOUT,7120) dt
 7120     format(10x,'TIME IN BACK SUBSTITUTION:',f9.5)
      endif
c
c normal exit
c------------
c
      return
c
c error exit
c ==========
c
 9100 continue
c
      return
      end




