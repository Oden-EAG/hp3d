      subroutine zdirio (Unname, Commnd, Irec, Len, Sbuf, Jerr)
c
#include "syscomz.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION: This routine preforms direct access I/O for scratch tapes
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Unname - name of the unit to operate on
c               (e.g. for surfs: 'U', 'B', or 'L')
c I   Commnd - character variable that tells what i/o action to take
c               it may be equal to 'OPEN', 'CLOSE', 'READ', or 'WRITE'
c I   Irec   - logical record number to read or write (not used for open
c
c I   Len    - for write: the number of words to write
c O   Len      for read: the number of words read
c I   Len      for open: the maximum length in words for a record
c
c I   Sbuf   - for write: array for writing out to
c O   Sbuf   - for read:  array for reading into
c
c O   Jerr   - error flag set to 0 for no error
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
#include "surfscwr.blk"
c
      character*(*) Unname, Commnd
cwb >
cwb       dimension Sbuf(*), nbuf(20), lenf(20)
      dimension Sbuf(*), nbuf(9), lenf(9)
      dimension lbuf(9), irsave(9)
c
cwb       save nbuf,lenf
ccccc save nbuf, lenf, irsave
      common /czdirio1/  nbuf, lenf, irsave   
      common /czdirio2/ storage(100)


c  ...test......test......test......test......test......test......test...
c      allocatable sbuf1(:)
c      character(6) TAPES(40)
c      data TAPES /'TAPE41','TAPE42','TAPE43','TAPE44','TAPE45',
c     .            'TAPE46','TAPE47','TAPE48','TAPE49','TAPE50',
c     .            'TAPE51','TAPE52','TAPE53','TAPE54','TAPE55',
c     .            'TAPE56','TAPE57','TAPE58','TAPE59','TAPE60',
c     .            'TAPE61','TAPE62','TAPE63','TAPE64','TAPE65',
c     .            'TAPE66','TAPE67','TAPE68','TAPE69','TAPE70',
c     .            'TAPE71','TAPE72','TAPE73','TAPE74','TAPE75',
c     .            'TAPE76','TAPE77','TAPE78','TAPE79','TAPE80'/
c  ...test......test......test......test......test......test......test...

cwb <
c
c the following data statement is machine dependant
c ...................................................
c within zdirio1.blk is contained the data statement:
c
      data maxrec,lenwrd /10000, 16/
c for complex(8), lenwrd = 16
c where:
c   maxrec = the maximum record length for a direct access file
c     **Note: on most machines this is bytes
c             on the cdc this is 8 byte words
c             on the ardent this is 4 byte words
c
c   lenwrd = the number of machine words in a floating point word
c     **Note: this actually indicates whether recl is in words or bytes
c            vax,ibm,apollo,etc (most other machines) = 8
c            cray,cdc = 1
c            ardent = 2
c
cccccccc       include '../newcomm/zdirio1.blk'
c ..................................................
c set the logical unit numbers used herein
c =========================================
cwb >
      data lbuf / 19,18,20, 21,22,23, 24,25,26 /
cwb       data lubufu, lubufb, lubufl / 19, 18, 20 /
cwb <

      if(IDUMPWR.eq.1)then
        write(nfsout,*)'UNEXPECTED CALL TO ZDIRIO - STOP'
        stop 'zdirio 1'
      endif

      irecp = 0
c
cwb >
c No resolution
c --------------
      ioffst = 0
c
c Unsymmetric w/ resolution capability
c ----------------------------------
      if (ISYM .eq. 3) ioffst = 3
c
c Symmetric w/ resolution capability
c ----------------------------------
      if (ISYM .eq. 4) ioffst = 6
c
c calculate the unit number from the unit name
c==============================================
      if (Unname .eq. 'U' ) then
         iunit = ioffst + 1
c
      elseif (Unname .eq. 'B' ) then
         iunit = ioffst + 2
c
      elseif (Unname .eq. 'L' ) then
         iunit = ioffst + 3
c
c unknown file
c
      else
         go to 9999
      endif
c
      ntape = lbuf(iunit)
c
cwb        if (Unname.eq.'U') then
cwb          iunit = 1
cwb          ntape = lubufu
cwb       elseif (Unname.eq.'B') then
cwb          iunit = 2
cwb          ntape = lubufb
cwb       elseif (Unname.eq.'L') then
cwb          iunit = 3
cwb          ntape = lubufl
cwb c
cwb c unknown file
cwb c
cwb        else
cwb          go to 9999
cwb        endif
cwb <
c-----------------------------------------------------------------------
c check for OPEN command
c ==========****========
c
      if (Commnd .eq. 'OPEN') then
c
c calculate the number of physical words from the logical words
c
         lenop = (Len+1)*lenwrd
c
c save the buffer size
c
         if (lenop.ge.maxrec) then
            nbuf(iunit) = 1 + lenop/maxrec
            lenop = maxrec
            lenf(iunit) = (maxrec/lenwrd) - 1
         else
            nbuf(iunit) = 1
            lenf(iunit) = Len
         endif
cwb >
         irsave(iunit) = 0
cwb <
c open the file
c
         close (ntape,err=5)
c
    5    continue
c
c the following open statement is machine dependant
c .................................................
c within zdirio2.blk is contained the apollo open statement:
c
         open (ntape,access='direct',form='unformatted',
     .         status='unknown',recl=lenop,err=9000)
c
ccccccccc       include '../newcomm/zdirio2.blk'
c..................................................
c
c normal return
c
         go to 1111
c-----------------------------------------------------------------------
c check for WRITE command
c ==========*****========
c
      elseif (Commnd.eq.'WRITE') then
cwb >
c Store the number of buffer records written
c   This is for resolution and setting the proper value of IFU,etc.
c
         irsave(iunit) = max(irsave(iunit),Irec)
cwb <
c calculate the first buffer to write
c
         irecp = 1 + (Irec-1)*nbuf(iunit)
         lenw = min0(lenf(iunit),Len)
c
c write the first buffer
c
         slen = Len
         write(ntape,rec=irecp,err=9200,iostat=iostat)
     .                            slen,(Sbuf(i),i=1,lenw)
         storage(ntape) = slen
c        ---------------------------------------------
c normal return
c
         if (Len .eq. lenw) go to 1111
c
c write out the rest of the records
c
         iend = lenw
c
 10      continue
c*****************
         irecp = irecp+1
         inow = iend+1
         iend = inow+lenf(iunit)
         if (iend .gt. Len) iend = Len
c
         write (ntape,rec=irecp,err=9200,iostat=iostat)
     .                               (Sbuf(i),i=inow,iend)
c ------------------------------------------------------
c normal return
c
         if (Len .eq. iend) go to 1111
c
         go to 10
c        --------
c-----------------------------------------------------------------------
c check for READ command
c===========****=========
c
      elseif (Commnd .eq. 'READ') then
cwb >
c Pick up the number of buffer records from the previous solution
c   This is for resolution and setting the proper value of IFU,etc.
c
         irecsv = irsave(iunit)
         if (Irec .lt. 0) Irec = irecsv
cwb <
c calculate the first buffer to read
c
         irecp = 1 + (Irec-1)*nbuf(iunit)
c
c read the buffer length
c
         slen=storage(ntape)
         read (ntape, rec=irecp, err=9300) slen
         Len = slen
         lenr = min0(lenf(iunit),Len)
c
c read the first buffer
c
         read (ntape, rec=irecp, err=9300) slen,(Sbuf(i),i=1,lenr)
c        ---------------------------------------------------------
c normal exit
c
         if (Len .eq. lenr) go to 1111
c
c read in the rest of the records
c
         iend = lenr
c
 20      continue
c****************
         irecp = irecp+1
         inow = iend+1
         iend = inow+lenf(iunit)
         if (iend .gt. Len) iend = Len
c
         read (ntape, rec=irecp, err=9300) (Sbuf(i),i=inow,iend)
c        --------------------------------------------------
c normal exit
c
         if (Len .eq. iend) go to 1111
c
         go to 20
c        --------
c-----------------------------------------------------------------------
c check for CLOSE
c ==========*****
c
      elseif (Commnd .eq. 'CLOSE') then
         close (ntape)
c        -------------
cwb >
         irsave(iunit) = 0
         nbuf(iunit) = 0
         lenf(iunit) = 0
cwb <
         go to 1111
      endif
c-----------------------------------------------------------------------
c unknown command
c ===============
c
      go to 9999
c
c NORMAL EXIT
c ===========
c
 1111 Jerr = 0


cwr07.12.00
c  ...test......test......test......test......test......test......test...
c  ...test......test......test......test......test......test......test...
c************************************************************************
c      if(commnd.eq.'WRITE')then
c      if(irec.gt.40)then
c        write(*,*)'ZDIRIO: IREC TOO BIG =',irec
c        stop
c      endif
c      write(*,*)'DUMPING...irec=',irec
c----------------------------
c        open(40,file=TAPES(irec))
c        write(40,*)(Sbuf(ii),ii=1,len)
c        close(40)
c      elseif(commnd.eq.'READ')then
c----------------------------
c        allocate( sbuf1(len), STAT = i1)
c        if(i1.ne.0)then
c          write(*,*)'ZDIRIO: did not allocate sbuf1'
c          stop
c        endif
c        open(40,file=TAPES(irec))
c        read(40,*)(Sbuf1(ii),ii=1,len)
c        close(40)
c    ...compare sbuf,sbuf1
c        do ii=1,len
c           if(abs(sbuf(ii)-sbuf1(ii)).gt.1.d-10)then
c             write(*,*)'ZDIRIO: sbuf1-sbuf=',abs(sbuf(ii)-sbuf1(ii))
c             write(*,*)'IREC,II=',irec,ii
c             stop
c           endif
c           enddo
c        deallocate(sbuf1,STAT = i1)
c      endif
c************************************************************************
c  ...test......test......test......test......test......test......test...
c  ...test......test......test......test......test......test......test...


c
c debug print
c
      if (IPFSZD .eq. 1) then
         write(NFSOUT,7701) Unname,Commnd,ntape,Irec,
     .                       Len,irecp,nbuf(iunit),lenf(iunit)
 7701     format(1x,'ZDIRIO: Unname,Commnd,ntape,Irec,',
     .           'Len,irecp,nbuf(iunit),lenf(iunit)',/,10x,a1,1x,a8,6i8)
      endif
      return
c
c-----------------------------------------------------------------------
c ERROR EXITS
c ===========
c open error
c
 9000 continue
      Jerr = 1
      write(NFSOUT,*) 'ERROR IN OPENING DIRECT ACCESS FILE'
      write(NFSOUT,7020) Unname,ntape
 7020  format(2(/),5x,'TAPE NAME',a8,2x,'UNIT - ',i3)
      return
c
c read error
c
 9200 continue
      Jerr = 2
      write(NFSOUT,*) 'ERROR IN WRITING TO DIRECT ACCESS FILE'
      write(NFSOUT,7040) Unname,ntape,Irec,Len,iostat
 7040  format(2(/),5x,
     .      '                BUFFER - ',a8,2x,'UNIT - ',i3,
     . /,5x,'                RECORD - ',i5,2x,'LENGTH - ',i7
     . /,5x,'                IOSTAT - ',i32)
      return
c
c write error
c
 9300 continue
      Jerr = 3
      write(NFSOUT,*) 'ERROR IN READING FROM DIRECT ACCESS FILE'
      write(NFSOUT,7060) Unname,ntape,Irec,Len
 7060 format(2(/),5x,
     .      '                BUFFER - ',a1,2x,'UNIT - ',i3,
     . /,5x,'                RECORD - ',i5,2x,'LENGTH - ',i7)
      return
c
c unknown command
c
 9999 Jerr = 5
      write(NFSOUT,*) 'IO PROBLEMS: UNKNOWN COMMAND'
      write(NFSOUT,7000) Unname,Commnd,Irec,Len,ntape
 7000  format(2(/),5x,
     . 'UNNAME,COMMND,IREC,LEN,NTAPE ',a1,2x,a8,2x,a1,a8,3i6)
      return
c
      end
