      subroutine read_sets
c
      use data_structure3D
      use control
c
#include "syscom.blk"
#include "cinout.blk"
#include "cnelcon.blk"
#include "csets.blk"
c
c  ...to read element sets 
      dimension istep(10),nrelrnd(10)
c

c  ...reading sets for graphical processing
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      NRSETS = 0


      read(LIN,*)NRSETS

      write(*,*) 'main: READING SETS, NRSETS = ',NRSETS
c
      if(NRSETS.gt.0)then
      ic=1
      do i=1,NRSETS

        NRELSET(2,i) = ic
        NRELSET(1,i) = 0
c
c  .....select case for the input file type
        read(LIN,*) iround,(nrelrnd(k),istep(k),k=1,iround)
c
c  .....loop through rounds:
        do k=1,iround
c
          if (istep(k).eq.0) then
c             
c  .........all elements are read:
              read(LIN,*)(NELISET(ic-1+j),j=1,nrelrnd(k))
c             
          else
             read(LIN,*)iaux

c  .........only first element is read, other found
            do j=1,nrelrnd(k)
              NELISET(ic-1+j) = iaux + (j-1)*istep(k)
            enddo
          endif
c
          ic = ic + nrelrnd(k)
          NRELSET(1,i) = NRELSET(1,i) + nrelrnd(k)
c
c  .....end of loop through rounds
        enddo

c      write(*,*)'read_sets, elset=',i
c      write(*,*)'nrelset(1:2)=',NRELSET(1:2,i) ; call pause
c      do ii=1,NRELSET(1,i)
c         write(*,*)ii,NELISET(NRELSET(2,i)+ii-1)
c      enddo
c      call pause
c
c  ...end of loop through sets
      enddo
      endif
c
      return
      end
