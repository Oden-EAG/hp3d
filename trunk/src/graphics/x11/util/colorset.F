c----------------------------------------------------------------------
c
c   routine name       - colorset
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Mar 02
c
c   purpose            - routine takes care of colors by creating
c                        a new colormap (when necessary) and next
c                        selecting entries for matrix NPCOL in order :
c                          NPCOL(1)-white,NPCOL(2)-black,
c                          NPCOL(3-10)-low-order-to-high-order-p-approx
c                          NPCOL(11-200) - spectrum of colors
c
c   usage              - call colorset(Iwind)
c
c   arguments :
c     in:
c               Iwind  = 1  white/black monitor
c                      = 2,3 gray scale monitor
c                      = 4 color monitor with 8 colors
c                      = 5 color monitors with 256 colors
c     out:
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine colorset(Iwind)
c
      use graphmod
#include "syscom.blk"
c
c.....colors to interpolate between
      parameter(nscales=5)
      dimension iram(4,nscales)
c
c  ...nscales=6
ccc      data iram /  0,  0,  0,255,   200,  0,255,255,
ccc     .           350,  0,255,  0,   500,255,255,  0,
ccc     .           800,255,  0,  0,  1000,255,  0,255/
c
c  ...from Tadeusz Liszka...
c  ...nscales=4
ccc      data iram /  0,  0,  0,255,  400,  0,255,  0,
ccc     .           700,255,255,  0, 1000,255,  0,  0/
c
c  ...nscales=5
      data iram /  0,  0,  0,255,  330,  0,255,255,
     .           500,  0,255,  0,  670,255,255,  0,
     .          1000,255,  0,  0/
c
c  ...nscales=9
ccc      data iram /  0,  0,  0,255,  100,  0,127,255,
ccc     .           250,  0,255,255,  450,  0,255,127,
ccc     .           500,  0,255,  0,  550,127,255,  0,
ccc     .           750,255,255,  0,  900,255,127,  0,
ccc     .          1000,255,  0,  0/
c
      iprint=0
c
c  ...first clear array of colors
      NPCOL(1:NR_COLORS) = 0
c
c  ...set white
      call xgsetcm(0,255,255,255)
      call   afill(1,255,255,255)
      NPCOL(1)=0
c
c  ...set black
      call xgsetcm(1,  0,  0,  0)
      call   afill(2,  0,  0,  0)
      NPCOL(2)=1
c
c  ...when on B&W machine
      if (Iwind.eq.1) return
c
c  ...when on gray scale
      if ((Iwind.eq.2).or.(Iwind.eq.3)) then
        ncol = 8
        ndcol = 255 / (ncol+1)
        do i = 1,ncol
          irgb = 255-ndcol*i
          call xgsetcm(i+1,irgb,irgb,irgb)
          NPCOL(i+2)=i+1
        enddo
        if (Iwind.eq.2) then
          do j=1,190
            NPCOL(j+10)=j/24+2
          enddo
c
c  .....when more colors possible on grayscale monitor
        elseif (Iwind.eq.3) then
          ncol = 190
          do i = 1,ncol
            irgb = 235-i
            call xgsetcm(i+9,irgb,irgb,irgb)
            NPCOL(i+10)=i+9
          enddo
        endif
c
c  ...when on color machine
      elseif ((Iwind.eq.4).or.(Iwind.eq.5)) then
c
c  .....define the colors for the p-scale...
c
c  .....p=1
        ir=0; ig=155; ib=255
        call xgsetcm(2,ir,ig,ib)
        call   afill(3,ir,ig,ib)
        NPCOL(3)=2
c
c  .....p=2
        ir=0; ig=255; ib=255
        call xgsetcm(3,ir,ig,ib)
        call   afill(4,ir,ig,ib)
        NPCOL(4)=3
c
c  .....p=3
        ir=0; ig=255; ib=0
        call xgsetcm(4,ir,ig,ib)
        call   afill(5,ir,ig,ib)
        NPCOL(5)=4
c
c  .....p=4
        ir=255; ig=255; ib=0
        call xgsetcm(5,ir,ig,ib)
        call   afill(6,ir,ig,ib)
        NPCOL(6)=5
c
c  .....p=5
        ir=255; ig=155; ib=0
        call xgsetcm(6,ir,ig,ib)
        call   afill(7,ir,ig,ib)
        NPCOL(7)=6
c
c  .....p=6
        ir=255; ig=25; ib=0
        call xgsetcm(7,ir,ig,ib)
        call   afill(8,ir,ig,ib)
        NPCOL(8)=7
c
c  .....p=7
        ir=255; ig=0; ib=155
        call xgsetcm(8,ir,ig,ib)
        call   afill(9,ir,ig,ib)
        NPCOL(9)=8
c
c  .....p=8
        ir=255; ig=0; ib=255
        call xgsetcm(9,ir,ig,ib)
        call  afill(10,ir,ig,ib)
        NPCOL(10)=9
c
c  .....define the colors for contour plots
        if (Iwind.eq.4) then
c
c  .......use the colors already defined
          do i=0,NR_COLORS-11
            NPCOL(11+i)=float(i)/float(NR_COLORS-10)*8 + 2
          enddo
        elseif (Iwind.eq.5) then
c
c  .......interpolate between the colors on 'iram' scale
          nrcol = 0
          nrrow = 1
c
c  .......loop through the colors to be generated
          do i=0,NR_COLORS-11
            if (nrcol.ge.iram(1,nrrow+1)) nrrow = nrrow + 1
            alpha = float(nrcol-iram(1,nrrow))/
     .                (iram(1,nrrow+1)-iram(1,nrrow))
            ir = alpha*iram(2,nrrow+1)+(1-alpha)*iram(2,nrrow)
            ig = alpha*iram(3,nrrow+1)+(1-alpha)*iram(3,nrrow)
            ib = alpha*iram(4,nrrow+1)+(1-alpha)*iram(4,nrrow)
c
            if (iprint.eq.1) then
              write(*,7003) nrcol,nrrow,alpha
 7003         format('nrcol,nrrow,alpha = ',2i4,f8.3)
              write(*,7001) i,ir,ig,ib
 7001         format('colorset: i = ',i4,' ir,ig,ib = ',3i4)
            endif
c
            call xgsetcm(10+i,ir,ig,ib)
            call   afill(11+i,ir,ig,ib)
            nrcol = nrcol + 1000/(NR_COLORS-10)
            NPCOL(11+i)=10+i
          enddo
c
        endif
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'COLORSET: NPCOL = '
        write(*,7002) NPCOL
 7002   format(10i4)
        call pause
      endif
c
      return
      end
c----------------------------------------------------------------------
c
c   routine name       - afill
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Feb 98
c
c   purpose            - an axiliary routine for colorset
c
c----------------------------------------------------------------------
c
      subroutine afill(index1,ir1,ig1,ib1)
c
      use graphmod
c
      IRED(index1)=ir1
      IGREEN(index1)=ig1
      IBLUE(index1)=ib1
c
      return
      end


