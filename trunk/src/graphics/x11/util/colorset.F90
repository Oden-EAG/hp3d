#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - colorset
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine takes care of colors by creating
!                        a new colormap (when necessary) and next
!                        selecting entries for matrix NPCOL in order :
!                          NPCOL(1)-white,NPCOL(2)-black,
!                          NPCOL(3-10)-low-order-to-high-order-p-approx
!                          NPCOL(11-200) - spectrum of colors
!
!   arguments :
!     in:
!               Iwind  = 1  white/black monitor
!                      = 2,3 gray scale monitor
!                      = 4 color monitor with 8 colors
!                      = 5 color monitors with 256 colors
!
!----------------------------------------------------------------------
      subroutine colorset(Iwind)
!
      use graphmod
!
      implicit none
!
      integer, intent(in) :: Iwind
!
      real(8) :: alpha
      integer :: i,ib,ig,ir,irgb,j
      integer :: ncol,ndcol,nrcol,nrrow
!
!.....colors to interpolate between
      integer, parameter :: nscales = 5
      integer :: iram(4,nscales)
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!  ...nscales=6
!!!      data iram /  0,  0,  0,255,   200,  0,255,255, &
!!!                 350,  0,255,  0,   500,255,255,  0, &
!!!                 800,255,  0,  0,  1000,255,  0,255/
!
!  ...from Tadeusz Liszka...
!  ...nscales=4
!!!      data iram /  0,  0,  0,255,  400,  0,255,  0, &
!!!                 700,255,255,  0, 1000,255,  0,  0/
!
!  ...nscales=5
      iram = reshape( (/  0,  0,  0,255,  330,  0,255,255, &
                        500,  0,255,  0,  670,255,255,  0, &
                       1000,255,  0,  0 /) , (/4,5/) )
!
!  ...nscales=9
!!!      data iram /  0,  0,  0,255,  100,  0,127,255, &
!!!                 250,  0,255,255,  450,  0,255,127, &
!!!                 500,  0,255,  0,  550,127,255,  0, &
!!!                 750,255,255,  0,  900,255,127,  0, &
!!!                1000,255,  0,  0/
!
!
!  ...first clear array of colors
      NPCOL(1:NR_COLORS) = 0
!
!  ...set white
      call xgsetcm(0,255,255,255)
      call   afill(1,255,255,255)
      NPCOL(1)=0
!
!  ...set black
      call xgsetcm(1,  0,  0,  0)
      call   afill(2,  0,  0,  0)
      NPCOL(2)=1
!
!  ...when on B&W machine
      if (Iwind.eq.1) return
!
!  ...when on gray scale
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
!
!  .....when more colors possible on grayscale monitor
        elseif (Iwind.eq.3) then
          ncol = 190
          do i = 1,ncol
            irgb = 235-i
            call xgsetcm(i+9,irgb,irgb,irgb)
            NPCOL(i+10)=i+9
          enddo
        endif
!
!  ...when on color machine
      elseif ((Iwind.eq.4).or.(Iwind.eq.5)) then
!
!  .....define the colors for the p-scale...
!
!  .....p=1
        ir=0; ig=155; ib=255
        call xgsetcm(2,ir,ig,ib)
        call   afill(3,ir,ig,ib)
        NPCOL(3)=2
!
!  .....p=2
        ir=0; ig=255; ib=255
        call xgsetcm(3,ir,ig,ib)
        call   afill(4,ir,ig,ib)
        NPCOL(4)=3
!
!  .....p=3
        ir=0; ig=255; ib=0
        call xgsetcm(4,ir,ig,ib)
        call   afill(5,ir,ig,ib)
        NPCOL(5)=4
!
!  .....p=4
        ir=255; ig=255; ib=0
        call xgsetcm(5,ir,ig,ib)
        call   afill(6,ir,ig,ib)
        NPCOL(6)=5
!
!  .....p=5
        ir=255; ig=155; ib=0
        call xgsetcm(6,ir,ig,ib)
        call   afill(7,ir,ig,ib)
        NPCOL(7)=6
!
!  .....p=6
        ir=255; ig=25; ib=0
        call xgsetcm(7,ir,ig,ib)
        call   afill(8,ir,ig,ib)
        NPCOL(8)=7
!
!  .....p=7
        ir=255; ig=0; ib=155
        call xgsetcm(8,ir,ig,ib)
        call   afill(9,ir,ig,ib)
        NPCOL(9)=8
!
!  .....p=8
        ir=255; ig=0; ib=255
        call xgsetcm(9,ir,ig,ib)
        call  afill(10,ir,ig,ib)
        NPCOL(10)=9
!
!  .....define the colors for contour plots
        if (Iwind.eq.4) then
!
!  .......use the colors already defined
          do i=0,NR_COLORS-11
            NPCOL(11+i)=int(float(i)/float(NR_COLORS-10)*8 + 2)
          enddo
        elseif (Iwind.eq.5) then
!
!  .......interpolate between the colors on 'iram' scale
          nrcol = 0
          nrrow = 1
!
!  .......loop through the colors to be generated
          do i=0,NR_COLORS-11
            if (nrcol.ge.iram(1,nrrow+1)) nrrow = nrrow + 1
            alpha = float(nrcol-iram(1,nrrow))/ &
                      (iram(1,nrrow+1)-iram(1,nrrow))
            ir = int(alpha*iram(2,nrrow+1)+(1-alpha)*iram(2,nrrow))
            ig = int(alpha*iram(3,nrrow+1)+(1-alpha)*iram(3,nrrow))
            ib = int(alpha*iram(4,nrrow+1)+(1-alpha)*iram(4,nrrow))
!
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,7003) nrcol,nrrow,alpha
 7003         format('nrcol,nrrow,alpha = ',2i4,f8.3)
              write(*,7001) i,ir,ig,ib
 7001         format('colorset: i = ',i4,' ir,ig,ib = ',3i4)
            endif
#endif
!
            call xgsetcm(10+i,ir,ig,ib)
            call   afill(11+i,ir,ig,ib)
            nrcol = nrcol + int(1000.e0/(NR_COLORS-10))
            NPCOL(11+i)=10+i
          enddo
!
        endif
      endif
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'COLORSET: NPCOL = '
        write(*,7002) NPCOL
 7002   format(10i4)
        call pause
      endif
#endif
!
      end subroutine colorset
!
!
!----------------------------------------------------------------------
!
!   routine name       - afill
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - an auxiliary routine for colorset
!
!----------------------------------------------------------------------
!
      subroutine afill(index1,ir1,ig1,ib1)
!
      use graphmod
!
      implicit none
!
      integer, intent(in) :: index1,ir1,ig1,ib1
!
      IRED  (index1)=ir1
      IGREEN(index1)=ig1
      IBLUE (index1)=ib1
!
!
      end subroutine afill

#endif
