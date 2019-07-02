      subroutine zdirio (Unname, Commnd, Irec, Len, Sbuf, Jerr)
c
#include "syscom.blk"
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
c
      character*(*) Unname, Commnd
cwb >
cwb       dimension Sbuf(*), nbuf(20), lenf(20)
      dimension Sbuf(*), nbuf(9), lenf(9)
      dimension lbuf(9), irsave(9)
c



cwb       save nbuf,lenf
      save nbuf, lenf, irsave
cwb <
c

c  ...new !!!
c     ~~~~~~~
c  ...irecpmax * maxrec = max files size in bytes
c  ...   10^5  * 10^4   = 10^9 = 1 Gbyte
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      save numfile, ntape0, lenop19
      data irecpmax /100000/
ccccc data irecpmax /5/


c the following data statement is machine dependant
c ...................................................
c within zdirio1.blk is contained the data statement:
c
      data maxrec,lenwrd /10000, 8/
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

c         write(*,*)'ZDIRIO: standard open, ntape=',ntape
c         call pause

         open (ntape,access='direct',form='unformatted',
     .                          recl=lenop,err=9000)
ccccc.         status='scratch',recl=lenop,err=9000)


      if(ntape.eq.19)then
        lenop19 = lenop
c        write(*,*)'saved lenop=',lenop19
c        call pause
      endif

c
ccccccccc       include '../newcomm/zdirio2.blk'
c..................................................
c
c normal return
c
         go to 1111
c-----------------------------------------------------------------------
c***********************************************************************
c***********************************************************************
c***********************************************************************
c-----------------------------------------------------------------------



      elseif (Commnd.eq.'WRITE') then

c     ...Store the number of buffer records written
c     ...This is for resolution and setting the proper value of IFU,etc.
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         irsave(iunit) = max(irsave(iunit),Irec)


c     ...calculate the first buffer to write
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         irecp = 1 + (Irec-1)*nbuf(iunit)
         lenw = min0(lenf(iunit),Len)


c     ...write the first buffer
c        ~~~~~~~~~~~~~~~~~~~~~~
         slen = Len

c---------------------------------------------------------------
c         write(ntape,rec=irecp,err=9200,iostat=iostat)
c     .                            slen,(Sbuf(i),i=1,lenw)
c---------------------------------------------------------------
         if(irecp.eq.1)then
           numfile = 1
           ntape0 = ntape
         endif

c     ...close old file, open new file
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if((irecp-1)/irecpmax*irecpmax.eq.(irecp-1) .and. irecp.ne.1)
     .     then
           numfile_old = (irecp-1)/irecpmax
           numfile = numfile_old+1

           ntape0_old =ntape + (numfile_old-1) + 11*min(numfile_old-1,1)
c           write(*,*)'1 ZDIRIO: new close, ntape0_old=',ntape0_old
           close(ntape0_old)

           ntape0 = ntape + (numfile-1) + 11*min(numfile-1,1)
c           write(*,*)'1 ZDIRIO: new open, ntape0=',ntape0
           open (ntape0,access='direct',form='unformatted',
     .                          recl=lenop19,err=9000)
c           call pause
         endif


         irecp0 = irecp - (numfile-1)*irecpmax

c         write(*,*)'1 ntape0=',ntape0
c         write(*,*)'1 irecp,irecp0,lenw=',irecp,irecp0,lenw
c         call pause
         write(ntape0,rec=irecp0,err=9200,iostat=iostat)
     .                            slen,(Sbuf(i),i=1,lenw)
c---------------------------------------------------------------

c     ...normal return
         if (Len .eq. lenw) go to 1111



c     ...write out the rest of the records
         iend = lenw
c
 10      continue
c*****************
         irecp = irecp+1
         inow = iend+1
         iend = inow+lenf(iunit)
         if (iend .gt. Len) iend = Len
c

c---------------------------------------------------------------
c         write (ntape,rec=irecp,err=9200,iostat=iostat)
c     .                               (Sbuf(i),i=inow,iend)
c---------------------------------------------------------------
c     ...close old file, open new file
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if((irecp-1)/irecpmax*irecpmax.eq.(irecp-1) .and. irecp.ne.1)
     .     then
           numfile_old = (irecp-1)/irecpmax
           numfile = numfile_old+1

           ntape0_old =ntape + (numfile_old-1) + 11*min(numfile_old-1,1)
c           write(*,*)'2 ZDIRIO: new close, ntape0_old=',ntape0_old
           close(ntape0_old)

           ntape0 = ntape + (numfile-1) + 11*min(numfile-1,1)
c           write(*,*)'2 ZDIRIO: new open, ntape0=',ntape0
           open (ntape0,access='direct',form='unformatted',
     .                          recl=lenop19,err=9000)
c           call pause

         endif


         irecp0 = irecp - (numfile-1)*irecpmax

c         write(*,*)'2 ntape0=',ntape0
c         write(*,*)'2 irecp,irecp0,lenw=',irecp,irecp0,lenw
c         call pause
         write(ntape0,rec=irecp0,err=9200,iostat=iostat)
     .                               (Sbuf(i),i=inow,iend)
c---------------------------------------------------------------



c     ...normal return
         if (Len .eq. iend) go to 1111


         go to 10
c        --------



c-----------------------------------------------------------------------
c check for READ command
c===========****=========
c
      elseif (Commnd .eq. 'READ') then

c     ...Pick up the number of buffer records from the previous solution
c     ...This is for resolution and setting the proper value of IFU,etc.
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         irecsv = irsave(iunit)
         if (Irec .lt. 0) Irec = irecsv

c     ...calculate the first buffer to read
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         irecp = 1 + (Irec-1)*nbuf(iunit)

         numfile_new = (irecp-1)/irecpmax + 1
c         write(*,*)'1 read, irecp=',irecp
c         write(*,*)'1 read, numfile,numfile_new=',numfile,numfile_new
c         call pause

c-------------------------------------------------------------
ccccc    read (ntape, rec=irecp, err=9300) slen,(Sbuf(i),i=1,lenr)
c         read (ntape, rec=irecp, err=9300) slen0
c         len0 = slen0
c         lenr0 =min0(lenf(iunit),len0)
c         read (ntape, rec=irecp, err=9300) slen0,(Sbuf0(i),i=1,lenr0)
c         write(*,*)'read sbuf0 ok...'
c         call pause

c-------------------------------------------------------------
c     ...close old file, open new file
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




         if(numfile.ne.numfile_new)then
           numfile_old = numfile
           numfile = numfile_new

           ntape0_old =ntape + (numfile_old-1) + 11*min(numfile_old-1,1)

c           write(*,*)'1 ZDIRIO: new close, ntape0=',ntape0_old
           close(ntape0_old)

           ntape0 = ntape + (numfile-1) + 11*min(numfile-1,1)

c           write(*,*)'1 ZDIRIO: new open, ntape0=',ntape0
           open (ntape0,access='direct',form='unformatted',
     .                          recl=lenop19,err=9000)
c           call pause

         endif


         irecp0 = irecp - (numfile-1)*irecpmax

c---------------------------------------------------------------

c     ...read the buffer length
         read (ntape0, rec=irecp0, err=9300) slen
         Len = slen
         lenr = min0(lenf(iunit),Len)


c     ...read the first buffer
c        ~~~~~~~~~~~~~~~~~~~~~
c         write(*,*)'1 read, ntape0,irecp,irecp0=',ntape0,irecp,irecp0
c         call pause
         read (ntape0, rec=irecp0, err=9300) slen,(Sbuf(i),i=1,lenr)




c     ...normal exit
         if (Len .eq. lenr) go to 1111


c     ...read in the rest of the records
         iend = lenr


 20      continue
c****************
         irecp = irecp+1
         inow = iend+1
         iend = inow+lenf(iunit)
         if (iend .gt. Len) iend = Len



c-----------------------------------------------------------------
ccccc    read (ntape, rec=irecp, err=9300) (Sbuf(i),i=inow,iend)
c         read (ntape, rec=irecp, err=9300) (Sbuf0(i),i=inow,iend)
c-----------------------------------------------------------------
c     ...close old file, open new file
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         numfile_new = (irecp-1)/irecpmax + 1

c         write(*,*)'2 read, irecp=',irecp
c         write(*,*)'2 read, numfile,numfile_new=',numfile,numfile_new
c         call pause

cwr09.05.03
ccccc    if(numfile.ne.mumfile_new)then
         if(numfile.ne.numfile_new)then

           numfile_old = numfile
           numfile = numfile_new

           ntape0_old =ntape + (numfile_old-1) + 11*min(numfile_old-1,1)
c           write(*,*)'2 ZDIRIO: new close, ntape0_old=',ntape0_old
           close(ntape0_old)

           ntape0 = ntape + (numfile-1) + 11*min(numfile-1,1)
c           write(*,*)'2 ZDIRIO: new open, ntape0=',ntape0
           open (ntape0,access='direct',form='unformatted',
     .                          recl=lenop19,err=9000)
c           call pause

         endif

         
         irecp0 = irecp - (numfile-1)*irecpmax

c         write(*,*)'2 read, ntape0, irecp,irecp0=',ntape0,irecp,irecp0 
         read (ntape0, rec=irecp0, err=9300) (Sbuf(i),i=inow,iend)
c---------------------------------------------------------------




c     ...normal exit
         if (Len .eq. iend) go to 1111
c
         go to 20
c        --------



c***********************************************************************
c***********************************************************************
c***********************************************************************
c-----------------------------------------------------------------------
c check for CLOSE
c ==========*****
c
      elseif (Commnd .eq. 'CLOSE') then

c         write(*,*)'ZDIRIO: standard close, ntape=',ntape
c         call pause
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
