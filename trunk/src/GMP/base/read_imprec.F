c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : this is a LEGACY routine. Do not use it, unless
c                    you revise it!      
c//////////////////////////////////////////////////////////////////////
c
c
c----------------------------------------------------------------------
c
c   routine name       - read_imprec
c
c----------------------------------------------------------------------
c
c   latest revision    - Jul 04
c
c   purpose            - routine reads in geometry data for an
c                        implicit rectangle; the routine is designed
c                        to cope with old input files where the
c                        connectivity information for the
c                        rectangle was not necessarily consistent
c                        with the information on surfaces
c                        forming the rectangle
c
c   arguments      
c     in:
c                Irec  - rectangle number
c
c---------------------------------------------------------------------
c
      subroutine read_imprec(Irec)
c
      use GMP
#include "syscom.blk"
#include "cinout.blk"
c
c  ...rectangle surface numbers read from the input file, 
      dimension nsurf(5)
c
c  ...surface numbers for an implicit curve
      dimension nsurfc(2)
c
c  ...index search results
      dimension nloc(2)
c
c  ...correspondence between numbering used for edge connectivties
c     and surface definition
      dimension nbij(4)
c
      i4(j) = j-(j-1)/4*4
c
c--------------------------------------------------------------------
c
      if (Irec.eq.16) then
        iprint=0
      else
      iprint=0
      endif
c
      if (iprint.eq.1) then
        write(*,7001) Irec
 7001   format('read_imprec: DEBUGGING FOR IMPLICIT RECTANGLE ',
     .          i4)
      endif
c
c  ...initiate 
      RECTANGLES(Irec)%Idata(1:5) = 0
      nbij(1:4) = 0
c
c  ...read surfaces defining the rectangle
      read(NIN,*) nsurf(1:5)
c
c  ...ignore the sign
      nsurf(1:5) = iabs(nsurf(1:5))
c
c  ...store the first surface on which the rectangle is located,
c     this must be a correct information
      RECTANGLES(Irec)%Idata(1) = nsurf(1)
c
c  ...loop through the curves defining the rectangle edges
      do i=1,4
c
c  .....use the connectivity info
        nc = iabs(RECTANGLES(Irec)%EdgeNo(i))
        if (CURVES(nc)%type.eq.'ImpCir') then
          nsurfc(1:2) = iabs(CURVES(nc)%Idata(1:2))
c
c  .......check if the surfaces defining the curve 
c         coincide with the surfaces defining the rectangle
          do j=1,2
            call locate(nsurfc(j),nsurf,5, nloc(j))
          enddo
          if (nloc(1).gt.1) then
            RECTANGLES(Irec)%Idata(1+i) = nsurfc(1)
            nbij(i) = nloc(1)-1
          endif
          if (nloc(2).gt.1) then
            RECTANGLES(Irec)%Idata(1+i) = nsurfc(2)
            nbij(i) = nloc(2)-1
          endif
        endif
      enddo
      if (iprint.eq.1) then
        write(*,7005) RECTANGLES(Irec)%Idata(1:5)
 7005   format('read_imprec: SURFACE DATA AFTER FIRST STEP = ',5i6)
        write(*,7006) nbij(1:4)
 7006   format('read_imprec: nbij = ',5i6)
      endif
c
c  ...check if at least one entry has been filled
      iflag=0
      do i=1,4
        if (nbij(i).ne.0) then
          iflag = 1
          exit
        endif
      enddo
      if (iprint.eq.1) then
        write(*,7007) i
 7007   format('read_imprec: HAVE FOUND SURFACE FOR EDGE NUMBER ',i1)
      endif
c
      select case(iflag)
c
c  ...at least one common surface has been found, attempt to
c     complete the info on surfaces
      case(1)
c
c  ...determine the compatibility of orientations
      nflag=0
      do j=1,3
        j1 = i4(i+j)
        j2 = i4(nbij(i)+j)
        if (RECTANGLES(Irec)%Idata(1+j1).eq.nsurf(1+j2)) nflag=1
      enddo
      do j=1,3
        j1 = i4(i+j)
        j2 = i4(nbij(i)+4-j)
        if (RECTANGLES(Irec)%Idata(1+j1).eq.nsurf(1+j2)) nflag=-1
      enddo
      if (iprint.eq.1) then
        write(*,7008) nflag
 7008   format('read_imprec: ORIENTATIONS COMPATIBILITY FLAG = ',i2) 
        call pause
      endif
      if (nflag.eq.0) then
        write(*,7002) Irec
 7002   format('input_imprec: INSUFFICIENT DATA TO DETERMINE ',
     .         'COMPATIBILITY OF ORIENTATIONS FOR RECTANGLE ', i4,
     .         ' ASSUMING COMPATIBILITY...')
        nflag=1
        call pause
      endif
c
c  ...complete the info on surfaces...
      i1 = i; i2 = nbij(i)
c
c  ...the orientations are compatible
      if (nflag.eq.1) then
        do j=0,3,1
          j1 = i4(i1+j); j2 = i4(i2+j)
          if (RECTANGLES(Irec)%Idata(1+j1).eq.0) 
     .        RECTANGLES(Irec)%Idata(1+j1) = nsurf(1+j2)
        enddo
c
c  ...the orientations are incompatible
      else
        do j=0,3,1
          j1 = i4(i1+j); j2 = i4(4+i2-j)
          if (RECTANGLES(Irec)%Idata(1+j1).eq.0) 
     .        RECTANGLES(Irec)%Idata(1+j1) = nsurf(1+j2)
        enddo
      endif
c
c  ...no common surfaces has been found
      case(0)
c
c  ...we have one more chance if all curves are implicit and
c     share the rectangle's first surface
c
c  ...loop through the curves defining the rectangle edges
      do i=1,4
c
c  .....use the connectivity info
        nc = iabs(RECTANGLES(Irec)%EdgeNo(i))
        if (CURVES(nc)%type.eq.'ImpCir') then
          nsurfc(1:2) = iabs(CURVES(nc)%Idata(1:2))
c
c  .......check if the surfaces defining the curve 
c         coincides= with the surfaces defining the rectangle
          do j=1,2
            call locate(nsurfc(j),nsurf,5, nloc(j))
          enddo
          if (nloc(1).eq.1) then
            RECTANGLES(Irec)%Idata(1+i) = nsurfc(2)
          endif
          if (nloc(2).eq.1) then
            RECTANGLES(Irec)%Idata(1+i) = nsurfc(1)
          endif
        endif
      enddo
c
c  ...cheack if all surface numbers have been filled
      nflag=1
      do i=1,4
        if (RECTANGLES(Irec)%Idata(1+i).eq.0) nflag=0
      enddo
c
      if (nflag.eq.0) then
        write(*,7003) Irec
 7003   format('input_imprec: INSUFFICIENT DATA TO DETERMINE ',
     .         'SURFACE DATA FOR FOR RECTANGLE ', i4,
     .         ' ASSUMING READ IN DATA')
        RECTANGLES(Irec)%Idata(2:5) = nsurf(2:5)
      endif
c
      end select
c
      if (iprint.eq.1) then
        write(*,7004) Irec, (RECTANGLES(Irec)%Idata(j),j=2,5),
     .                (nsurf(j),j=2,5)
 7004   format('read_imprec: Irec, RECONSTRUCTED AND ',
     .         'INPUTED SURFACES = ',i3,2x,4i4,2x,4i4)
        call pause
      endif
c
c
      end
