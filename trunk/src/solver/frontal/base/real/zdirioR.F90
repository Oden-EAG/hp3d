      
!
!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: This routine performs direct access I/O for scratch tapes
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   Unname - name of the unit to operate on
!               (e.g. for surfs: 'U', 'B', or 'L')
! I   Commnd - character variable that tells what i/o action to take
!               it may be equal to 'OPEN', 'CLOSE', 'READ', or 'WRITE'
! I   Irec   - logical record number to read or write (not used for open
!
! I   Len    - for write: the number of words to write
! O   Len      for read: the number of words read
! I   Len      for open: the maximum length in words for a record
!
! I   Sbuf   - for write: array for writing out to
! O   Sbuf   - for read:  array for reading into
!
! O   Jerr   - error flag set to 0 for no error
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
   subroutine zdirio (Unname, Commnd, Irec, Len, Sbuf, Jerr)
!
      use surfsc1
!
      implicit none
!
      character*(*) Unname, Commnd
!
      integer :: Irec, Len, Jerr
      real(8) :: Sbuf(*)
!
      integer :: nbuf(9), lenf(9)
      integer :: lbuf(9), irsave(9)
!
      integer :: slen
      integer :: i,iend,inow,ioffst,iostat,iunit
      integer :: irecp,irecp0,irecpmax,irecsv
      integer :: lenop,lenop19,lenr,lenw,lenwrd
      integer :: maxrec,ntape,ntape0,ntape0_old
      integer :: numfile,numfile_new,numfile_old
!
      save nbuf, lenf, irsave
!

!  ...new !!!
!     ~~~~~~~
!  ...irecpmax * maxrec = max files size in bytes
!  ...   10^5  * 10^4   = 10^9 = 1 Gbyte
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      save numfile, ntape0, lenop19
      data irecpmax /100000/
!!!!! data irecpmax /5/


! the following data statement is machine dependent
! ...................................................
!
      data maxrec,lenwrd /10000, 8/
! where:
!   maxrec = the maximum record length for a direct access file
!     **Note: on most machines this is bytes
!             on the cdc this is 8 byte words
!             on the ardent this is 4 byte words
!
!   lenwrd = the number of machine words in a floating point word
!     **Note: this actually indicates whether recl is in words or bytes
!            vax,ibm,apollo,etc (most other machines) = 8
!            cray,cdc = 1
!            ardent = 2
!
! ..................................................
! set the logical unit numbers used herein
! =========================================
      data lbuf / 19,18,20, 21,22,23, 24,25,26 /
      irecp = 0
!
! No resolution
! --------------
      ioffst = 0
!
! Unsymmetric w/ resolution capability
! ----------------------------------
      if (ISYM .eq. 3) ioffst = 3
!
! Symmetric w/ resolution capability
! ----------------------------------
      if (ISYM .eq. 4) ioffst = 6
!
! calculate the unit number from the unit name
!==============================================
      if (Unname .eq. 'U' ) then
         iunit = ioffst + 1
!
      elseif (Unname .eq. 'B' ) then
         iunit = ioffst + 2
!
      elseif (Unname .eq. 'L' ) then
         iunit = ioffst + 3
!
! unknown file
!
      else
         go to 9999
      endif
!
      ntape = lbuf(iunit)
!
!wb        if (Unname.eq.'U') then
!wb          iunit = 1
!wb          ntape = lubufu
!wb       elseif (Unname.eq.'B') then
!wb          iunit = 2
!wb          ntape = lubufb
!wb       elseif (Unname.eq.'L') then
!wb          iunit = 3
!wb          ntape = lubufl
!wb !
!wb ! unknown file
!wb !
!wb        else
!wb          go to 9999
!wb        endif
!wb <
!-----------------------------------------------------------------------
! check for OPEN command
! ==========****========
!
      if (Commnd .eq. 'OPEN') then
!
! calculate the number of physical words from the logical words
!
         lenop = (Len+1)*lenwrd
!
! save the buffer size
!
         if (lenop.ge.maxrec) then
            nbuf(iunit) = 1 + lenop/maxrec
            lenop = maxrec
            lenf(iunit) = (maxrec/lenwrd) - 1
         else
            nbuf(iunit) = 1
            lenf(iunit) = Len
         endif
!wb >
         irsave(iunit) = 0
!wb <
! open the file
!
         close (ntape,err=5)
!
    5    continue
!
! the following open statement is machine dependant
! .................................................
!
!         write(*,*)'ZDIRIO: standard open, ntape=',ntape
!         call pause

         open (ntape,access='direct',form='unformatted', &
                                recl=lenop,err=9000)
!
      if(ntape.eq.19)then
        lenop19 = lenop
!        write(*,*)'saved lenop=',lenop19
!        call pause
      endif
!
!..................................................
!
! normal return
!
         go to 1111
!-----------------------------------------------------------------------
!***********************************************************************
!***********************************************************************
!***********************************************************************
!-----------------------------------------------------------------------



      elseif (Commnd.eq.'WRITE') then

!      ..Store the number of buffer records written
!      ..This is for resolution and setting the proper value of IFU,etc.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         irsave(iunit) = max(irsave(iunit),Irec)


!      ..calculate the first buffer to write
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         irecp = 1 + (Irec-1)*nbuf(iunit)
         lenw = min0(lenf(iunit),Len)


!      ..write the first buffer
!        ~~~~~~~~~~~~~~~~~~~~~~
         slen = Len

!---------------------------------------------------------------
!         write(ntape,rec=irecp,err=9200,iostat=iostat) &
!                                  slen,(Sbuf(i),i=1,lenw)
!---------------------------------------------------------------
         if(irecp.eq.1)then
           numfile = 1
           ntape0 = ntape
         endif

!      ..close old file, open new file
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if((irecp-1)/irecpmax*irecpmax.eq.(irecp-1) .and. irecp.ne.1) &
           then
           numfile_old = (irecp-1)/irecpmax
           numfile = numfile_old+1

           ntape0_old =ntape + (numfile_old-1) + 11*min(numfile_old-1,1)
!           write(*,*)'1 ZDIRIO: new close, ntape0_old=',ntape0_old
           close(ntape0_old)

           ntape0 = ntape + (numfile-1) + 11*min(numfile-1,1)
!           write(*,*)'1 ZDIRIO: new open, ntape0=',ntape0
           open (ntape0,access='direct',form='unformatted', &
                                recl=lenop19,err=9000)
!           call pause
         endif


         irecp0 = irecp - (numfile-1)*irecpmax

!         write(*,*)'1 ntape0=',ntape0
!         write(*,*)'1 irecp,irecp0,lenw=',irecp,irecp0,lenw
!         call pause
         write(ntape0,rec=irecp0,err=9200,iostat=iostat) &
                                  slen,(Sbuf(i),i=1,lenw)
!---------------------------------------------------------------

!      ..normal return
         if (Len .eq. lenw) go to 1111



!      ..write out the rest of the records
         iend = lenw
!
 10      continue
!*****************
         irecp = irecp+1
         inow = iend+1
         iend = inow+lenf(iunit)
         if (iend .gt. Len) iend = Len
!

!---------------------------------------------------------------
!         write (ntape,rec=irecp,err=9200,iostat=iostat) &
!                                     (Sbuf(i),i=inow,iend)
!---------------------------------------------------------------
!      ..close old file, open new file
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if((irecp-1)/irecpmax*irecpmax.eq.(irecp-1) .and. irecp.ne.1) &
           then
           numfile_old = (irecp-1)/irecpmax
           numfile = numfile_old+1

           ntape0_old =ntape + (numfile_old-1) + 11*min(numfile_old-1,1)
!           write(*,*)'2 ZDIRIO: new close, ntape0_old=',ntape0_old
           close(ntape0_old)

           ntape0 = ntape + (numfile-1) + 11*min(numfile-1,1)
!           write(*,*)'2 ZDIRIO: new open, ntape0=',ntape0
           open (ntape0,access='direct',form='unformatted', &
                                recl=lenop19,err=9000)
!           call pause

         endif


         irecp0 = irecp - (numfile-1)*irecpmax

!         write(*,*)'2 ntape0=',ntape0
!         write(*,*)'2 irecp,irecp0,lenw=',irecp,irecp0,lenw
!         call pause
         write(ntape0,rec=irecp0,err=9200,iostat=iostat) &
                                     (Sbuf(i),i=inow,iend)
!---------------------------------------------------------------



!      ..normal return
         if (Len .eq. iend) go to 1111


         go to 10
!        --------



!-----------------------------------------------------------------------
! check for READ command
!===========****=========
!
      elseif (Commnd .eq. 'READ') then

!      ..Pick up the number of buffer records from the previous solution
!      ..This is for resolution and setting the proper value of IFU,etc.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         irecsv = irsave(iunit)
         if (Irec .lt. 0) Irec = irecsv

!      ..calculate the first buffer to read
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         irecp = 1 + (Irec-1)*nbuf(iunit)

         numfile_new = (irecp-1)/irecpmax + 1
!         write(*,*)'1 read, irecp=',irecp
!         write(*,*)'1 read, numfile,numfile_new=',numfile,numfile_new
!         call pause

!-------------------------------------------------------------
!!!!!    read (ntape, rec=irecp, err=9300) slen,(Sbuf(i),i=1,lenr)
!         read (ntape, rec=irecp, err=9300) slen0
!         len0 = slen0
!         lenr0 =min0(lenf(iunit),len0)
!         read (ntape, rec=irecp, err=9300) slen0,(Sbuf0(i),i=1,lenr0)
!         write(*,*)'read sbuf0 ok...'
!         call pause

!-------------------------------------------------------------
!      ..close old file, open new file
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




         if(numfile.ne.numfile_new)then
           numfile_old = numfile
           numfile = numfile_new

           ntape0_old =ntape + (numfile_old-1) + 11*min(numfile_old-1,1)

!           write(*,*)'1 ZDIRIO: new close, ntape0=',ntape0_old
           close(ntape0_old)

           ntape0 = ntape + (numfile-1) + 11*min(numfile-1,1)

!           write(*,*)'1 ZDIRIO: new open, ntape0=',ntape0
           open (ntape0,access='direct',form='unformatted', &
                                recl=lenop19,err=9000)
!           call pause

         endif


         irecp0 = irecp - (numfile-1)*irecpmax

!---------------------------------------------------------------

!      ..read the buffer length
         read (ntape0, rec=irecp0, err=9300) slen
         Len = slen
         lenr = min0(lenf(iunit),Len)


!      ..read the first buffer
!        ~~~~~~~~~~~~~~~~~~~~~
!         write(*,*)'1 read, ntape0,irecp,irecp0=',ntape0,irecp,irecp0
!         call pause
         read (ntape0, rec=irecp0, err=9300) slen,(Sbuf(i),i=1,lenr)




!      ..normal exit
         if (Len .eq. lenr) go to 1111


!      ..read in the rest of the records
         iend = lenr


 20      continue
!****************
         irecp = irecp+1
         inow = iend+1
         iend = inow+lenf(iunit)
         if (iend .gt. Len) iend = Len



!-----------------------------------------------------------------
!!!!!    read (ntape, rec=irecp, err=9300) (Sbuf(i),i=inow,iend)
!         read (ntape, rec=irecp, err=9300) (Sbuf0(i),i=inow,iend)
!-----------------------------------------------------------------
!      ..close old file, open new file
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         numfile_new = (irecp-1)/irecpmax + 1

!         write(*,*)'2 read, irecp=',irecp
!         write(*,*)'2 read, numfile,numfile_new=',numfile,numfile_new
!         call pause

!wr09.05.03
!!!!!    if(numfile.ne.mumfile_new)then
         if(numfile.ne.numfile_new)then

           numfile_old = numfile
           numfile = numfile_new

           ntape0_old =ntape + (numfile_old-1) + 11*min(numfile_old-1,1)
!           write(*,*)'2 ZDIRIO: new close, ntape0_old=',ntape0_old
           close(ntape0_old)

           ntape0 = ntape + (numfile-1) + 11*min(numfile-1,1)
!           write(*,*)'2 ZDIRIO: new open, ntape0=',ntape0
           open (ntape0,access='direct',form='unformatted', &
                                recl=lenop19,err=9000)
!           call pause

         endif


         irecp0 = irecp - (numfile-1)*irecpmax

!         write(*,*)'2 read, ntape0, irecp,irecp0=',ntape0,irecp,irecp0
         read (ntape0, rec=irecp0, err=9300) (Sbuf(i),i=inow,iend)
!---------------------------------------------------------------




!      ..normal exit
         if (Len .eq. iend) go to 1111
!
         go to 20
!        --------



!***********************************************************************
!***********************************************************************
!***********************************************************************
!-----------------------------------------------------------------------
! check for CLOSE
! ==========*****
!
      elseif (Commnd .eq. 'CLOSE') then

!         write(*,*)'ZDIRIO: standard close, ntape=',ntape
!         call pause
         close (ntape)
!        -------------
!wb >
         irsave(iunit) = 0
         nbuf(iunit) = 0
         lenf(iunit) = 0
!wb <
         go to 1111
      endif
!-----------------------------------------------------------------------
! unknown command
! ===============
!
      go to 9999
!
! NORMAL EXIT
! ===========
!
 1111 Jerr = 0
!
! debug print
!
      if (IPFSZD .eq. 1) then
         write(NFSOUT,7701) Unname,Commnd,ntape,Irec, &
                             Len,irecp,nbuf(iunit),lenf(iunit)
 7701     format(1x,'ZDIRIO: Unname,Commnd,ntape,Irec,', &
                 'Len,irecp,nbuf(iunit),lenf(iunit)',/,10x,a1,1x,a8,6i8)
      endif
      return
!
!-----------------------------------------------------------------------
! ERROR EXITS
! ===========
! open error
!
 9000 continue
      Jerr = 1
      write(NFSOUT,*) 'ERROR IN OPENING DIRECT ACCESS FILE'
      write(NFSOUT,7020) Unname,ntape
 7020  format(2(/),5x,'TAPE NAME',a8,2x,'UNIT - ',i3)
      return
!
! read error
!
 9200 continue
      Jerr = 2
      write(NFSOUT,*) 'ERROR IN WRITING TO DIRECT ACCESS FILE'
      write(NFSOUT,7040) Unname,ntape,Irec,Len,iostat
 7040  format(2(/),5x, &
            '                BUFFER - ',a8,2x,'UNIT - ',i3, &
       /,5x,'                RECORD - ',i5,2x,'LENGTH - ',i7 &
       /,5x,'                IOSTAT - ',i32)
      return
!
! write error
!
 9300 continue
      Jerr = 3
      write(NFSOUT,*) 'ERROR IN READING FROM DIRECT ACCESS FILE'
      write(NFSOUT,7060) Unname,ntape,Irec,Len
 7060 format(2(/),5x, &
            '                BUFFER - ',a1,2x,'UNIT - ',i3, &
       /,5x,'                RECORD - ',i5,2x,'LENGTH - ',i7)
      return
!
! unknown command
!
 9999 Jerr = 5
      write(NFSOUT,*) 'IO PROBLEMS: UNKNOWN COMMAND'
      write(NFSOUT,7000) Unname,Commnd,Irec,Len,ntape
 7000  format(2(/),5x, &
       'UNNAME,COMMND,IREC,LEN,NTAPE ',a1,2x,a8,2x,a1,a8,3i6)
      return
!
!
   end subroutine zdirio
