      subroutine resol (Awrk)
c
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c initiate forward elimination of rhs only,  followed by backsubstitutio
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c W   Awrk  : workspace used by the solver (it is cut and passed below)
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
c
#include "syscom.blk"
c
      real(8) :: Awrk(*)
c
      call zecond(t0)
cwb >
      NNEGP = -1
      NPOSP = -1
c
      IFU = -1
      IFL = -1
cwb <
      IERR = 0
      IFB = 0
c
c check for invalid circumstances
c -------------------------------
      if (ISYM.eq.2 .or. ISYM.eq.1) then
         IERR=3
         write(NFSOUT,7020)
 7020    format(2(/),5x,
     .     'ERROR IN SURFS: ATTEMPTED RESOLUTION WITH',
     .     ' RESOLUTION INACTIVATED (IE: ISYM=2 OR 1)')
         return
      endif
c
      if (NRHS .eq. 0) then
         IERR=4
         write (NFSOUT,7040)
 7040    format(2(/), 5x,
     .     'ERROR IN SURFS: ATTEMPTED RESOLUTION WITH 0 RHS')
         return
      endif
c
c n = amount of space required to hold the first few temp arrays cut fro
c      the workspace
c
      n = NUMELM + MLDEST + 2*MDOF
c
c
cwb >
c  Pull up the old values of MBUF & MW to check against for resolution
c *** Note: we cannot let the size of MBUF shrink during resolution
c            because then the stuff that was buffered out will not fit b
c
      if (ISYM .eq. 3) then
         MBUF = MBUFSV(1,1)
         MW = MBUFSV(1,2)
      elseif (ISYM .eq. 4) then
         MBUF = MBUFSV(2,1)
         MW = MBUFSV(2,2)
      endif
c
c check against the current size of the workspace
c
      machk = MBUF + MW + n
      if (MA .lt. machk) then
         IERR=6
         write(NFSOUT,7027)
 7027    format(2(/),5x,
     .     'ERROR IN SURFS: ATTEMPTED RESOLUTION WITH',
     .     ' LESS WORKSPACE THAN WAS PREVIOUSLY ALLOCATED')
         return
      endif
c
cwb <
c
c
c set up storage cut lengths
c---------------------------
c assembly space cut length for elem rhs assembly
c
      MELEM = MDOF*NRHS*IASSEM
c
c cut length for the eleminated rhs equations
c
      MFWR = MFW*NRHS
c
c space which is available for buffering rhs equations
c
      MB = MW - MELEM - MFWR
c
c open direct access buffer files
c
      call zdirio ('B', 'OPEN', 0 , MB, Awrk, IERR)
c     -----------------------------------------
c error out
c
      IERR = 10*IERR
      if (IERR .ne. 0) go to 9100
c
c debug print
c
      if(IPRSTR.eq.1) then
         if(ISYM .eq. 4) then
c           write(NFSOUT,7060) NRHS
         else
c           write(NFSOUT,7080) NRHS
         endif
c        write(NFSOUT,7100) n,MELEM,MFWR,MB,MBUF,MA
c
 7060    format(2(/), 5x,'SYMMETRIC RESOLUTION WITH',i3,' RHS',/)
 7080    format(2(/), 5x,'UNSYMMETRIC RESOLUTION WITH',i3,' RHS',/)
 7100    format(     4x,'          OVER HEAD:',i7,/,
     .               4x,'            ELEMENT:',i7,/,
     .               4x,'              FRONT:',i7,/,
     .               4x,'         RHS BUFFER:',i7,/,
     .               4x,'         LHS BUFFER:',i7,/,
     .               4x,'      TOTAL STORAGE:',i7)
      endif
c
c not enuf space to buffer even 1 equation
c
      if (MB .lt. NRHS) then
         IERR = 5
         minn = max0(0,MB)
         write (NFSOUT,7120) minn
 7120    format(2(/),5x,'ERROR: TOO MANY RHS',
     .             /,12x,'MAXIMUM NUMBER OF RHS =',i2)
         go to 9100
      endif
c
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
      ial = ia1 + NUMELM
c
c set pointer to temp storage for dof destination vectors         (amdes
c
      iam = ial + MLDEST
c
c set pointer to temp storage for the transfer loc for active dof
c   in the front after eliminating                                (andes
c
      ian = iam + MDOF
c
c set pointer to temp storage for elem lhs + rhs, during assembly   (ele
c
      iae = ian + MDOF
c
c set pointer to temp storage for rhs equations for front           (frn
c
      iaf = iae + MELEM
c
c set pointer to temp storage for buffered rhs equations            (bbu
c
      iabr = iaf + MFWR
c
c set pointer to temp storage for equation buffering,
c   this is all remaining workspace                                  (bu
c
      iabf = iabr + MB
c
c
c forward elimination of all elements  (rhs only)
c ====================================
c
      call frwrs (Awrk(ia1), Awrk(ial), Awrk(iam), Awrk(ian),
     .            Awrk(iae), Awrk(iaf), Awrk(iabr),Awrk(iabf))
c     ---------------------------------------------
c
      call zecond(tf)
      dt = tf-t0
c
      if(IPRSTR .eq. 1) then
         write(NFSOUT,7140) dt
         if(IFB .gt. 0) write(NFSOUT,7160) IFB
 7140     format(10x,'TIME IN FORWARD ELIMINATION:',f9.3)
 7160     format(10x,'RHS BUFFER DUMPS:',i4)
      endif
c
      if (IERR .ne .0) go to 9100
c
c backward substitution of all elements
c ====================================
c
      call bckwrd ( Awrk(ia1), Awrk(ial), Awrk(iam), Awrk(ian),
     .              Awrk(iae) ,Awrk(iaf), Awrk(iabr),Awrk(iabf))
c      ----------------------------------------------
c
      call zecond(tb)
      dt = tb - t0
      if(IPRSTR .eq. 1) write(NFSOUT,7180) dt
 7180   format(10x,'TIME IN BACK SUBSTITUTION:',f9.3)
c
c normal exit
c------------
c
      return
c
c error exit
c ----------
 9100 continue
c
c
      end subroutine resol
