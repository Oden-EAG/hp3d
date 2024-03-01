!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION: This routine preforms direct access I/O for scratch tapes
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
   subroutine zdiriodum (Unname, Commnd, Irec, Len, Sbuf, Jerr)
!
      use czdirio
      use surfsc1
!
      implicit none
!
      character*(*) Unname, Commnd
!
      integer    :: Irec,Len,Jerr
      complex(8) :: Sbuf(*)
!
      integer    :: iend,inow,ioffst,iostat,irecp,irecsv,iunit
      integer    :: lenop,lenr,lenw,ntape
      complex(8) :: slen
      integer    :: lbuf(9)
!
! the following data statement is machine dependant
! ...................................................
!
      integer, parameter :: maxrec = 10000
      integer, parameter :: lenwrd = 16
! for complex(8), lenwrd = 16
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
         goto 9999
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
!wb          goto 9999
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
         open (ntape,access='direct',form='unformatted', &
               status='unknown',recl=lenop,err=9000)
!
!..................................................
!
! normal return
!
         goto 1111
!-----------------------------------------------------------------------
! check for WRITE command
! ==========*****========
!
      elseif (Commnd.eq.'WRITE') then
!wb >
! Store the number of buffer records written
!   This is for resolution and setting the proper value of IFU,etc.
!
         irsave(iunit) = max(irsave(iunit),Irec)
!wb <
! calculate the first buffer to write
!
         irecp = 1 + (Irec-1)*nbuf(iunit)
         lenw = min0(lenf(iunit),Len)
!
! write the first buffer
!
         slen = Len
!!!!!    write(ntape,rec=irecp,err=9200,iostat=iostat)
!!!!!.                            slen,(Sbuf(i),i=1,lenw)
         storage(ntape) = slen
!        ---------------------------------------------
! normal return
!
         if (Len .eq. lenw) goto 1111
!
! write out the rest of the records
!
         iend = lenw
!
 10      continue
!*****************
         irecp = irecp+1
         inow = iend+1
         iend = inow+lenf(iunit)
         if (iend .gt. Len) iend = Len
!
!!!!!    write (ntape,rec=irecp,err=9200,iostat=iostat)
!!!!!.                               (Sbuf(i),i=inow,iend)
! ------------------------------------------------------
! normal return
!
         if (Len .eq. iend) goto 1111
!
         goto 10
!        --------
!-----------------------------------------------------------------------
! check for READ command
!===========****=========
!
      elseif (Commnd .eq. 'READ') then
!wb >
! Pick up the number of buffer records from the previous solution
!   This is for resolution and setting the proper value of IFU,etc.
!
         irecsv = irsave(iunit)
         if (Irec .lt. 0) Irec = irecsv
!wb <
! calculate the first buffer to read
!
         irecp = 1 + (Irec-1)*nbuf(iunit)
!
! read the buffer length
!
         slen=storage(ntape)
!!!!!    read (ntape, rec=irecp, err=9300) slen
         Len = int(slen)
         lenr = min0(lenf(iunit),Len)
!
! read the first buffer
!
!!!!!    read (ntape, rec=irecp, err=9300) slen,(Sbuf(i),i=1,lenr)
!        ---------------------------------------------------------
! normal exit
!
         if (Len .eq. lenr) goto 1111
!
! read in the rest of the records
!
         iend = lenr
!
 20      continue
!****************
         irecp = irecp+1
         inow = iend+1
         iend = inow+lenf(iunit)
         if (iend .gt. Len) iend = Len
!
!!!!!    read (ntape, rec=irecp, err=9300) (Sbuf(i),i=inow,iend)
!        --------------------------------------------------
! normal exit
!
         if (Len .eq. iend) goto 1111
!
         goto 20
!        --------
!-----------------------------------------------------------------------
! check for CLOSE
! ==========*****
!
      elseif (Commnd .eq. 'CLOSE') then
         close (ntape)
!        -------------
!wb >
         irsave(iunit) = 0
         nbuf(iunit) = 0
         lenf(iunit) = 0
!wb <
         goto 1111
      endif
!-----------------------------------------------------------------------
! unknown command
! ===============
!
      goto 9999
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
!
!
   end subroutine zdiriodum
