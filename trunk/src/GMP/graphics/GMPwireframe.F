#if HP3D_USE_X11

c----------------------------------------------------------------------
c
c   routine name       - GMPwireframe
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine displays the wireframe for a 3D manifold
c
c
c   arguments :
c     in:
c               Iwind  - window number
c
c----------------------------------------------------------------------
c
      subroutine GMPwireframe(Iwind)
c
      use graphmod
      use GMP
c
      implicit none
c
      integer :: Iwind
c
c  ...list of segments
      real(8), allocatable ::  SEGMENTS(:,:,:)
      integer, allocatable ::  COLORN(:)
c
c  ...list of rectangles to display
      integer, parameter :: max_nr_rect=10
      integer :: nrect(max_nr_rect)
c
c  ...list of curves to display
      integer, parameter :: max_nr_curv=10
      integer :: ncurv(max_nr_curv)
c
c  ...point, derivatives
      real(8) :: x(3),dxdeta(3),xob(3),xyp(2,2)
c
c  ...bounds
      real(8) :: xmax(2),xmin(2)
      integer :: im(2)
c
      real(8) :: eta,psi,q,q1,rnsc,theta
      integer :: ic,ichoice,idec,idec_bundle,idec_curve,idec_edge
      integer :: idec_Hermite,idec_normal,idec_rectangle
      integer :: ii,ip,ir,j,k,lab,nc,ncol,np,npoint
      integer :: nr,nr_curv,nr_curves,nr_rect,nr_segments,nrsub1
c
      integer :: iprint
c
      real(8) :: bigp,bign,small,pi
      data bigp,bign,small /1.d30,-1.d30,1.d-14/
      data pi /3.14159265358979312d0/
c
c***********************************************************************
      iprint=0
c
c  ...set up default values for control parameters
      idec_Hermite=0
      idec_rectangle=0
      idec_curve=0
      idec_bundle=0
      idec_normal=0
      idec_edge=0
c
c  ...set initial scaling constants
      XCIM(1) = 0.d0
      XCIM(2) = 0.d0
      xmax(1:2) = bign
      xmin(1:2) = bigp
      XCWIN(1) = rmargin + xlength/2.d0
      XCWIN(2) = rmargin + ylength/2.d0
c
c  ...set up number of line segments for drawing a curve
      nrsub1 = 2
      NRSUB = 2**nrsub1
      DX = 1.d0/NRSUB
c
c
c  ...determine data for projection...
      theta = pi/4.d0
      psi = pi/4.d0
      RN(1) = cos(psi)*cos(theta)
      RN(2) = cos(psi)*sin(theta)
      RN(3) = sin(psi)
c
 10   continue
c
c  ...prepare data for transformation from physical coordinates to
c     the observer's system
      call cartobs
c
c  ...create the list of segments to display
      nr_segments = NRCURVE*NRSUB
      if (allocated(SEGMENTS)) deallocate(SEGMENTS)
      if (allocated(COLORN)) deallocate(COLORN)
      allocate(SEGMENTS(2,2,nr_segments))
      allocate(COLORN(nr_segments))
c
      ic=0; nr_curves=0
      do nc=1,NRCURVE
c
        if (iprint.eq.1) then
          write(*,7051) nc
 7051     format('GMPwireframe: nc = ',i6)
        endif
c
        ncol=NPCOL(2)
        call check_edge(nc, ii)
        if (ii.eq.2) ncol=NPCOL(8)
c
        if (idec_edge.eq.1) then
          call check_edge(nc, ii)
          if (ii.ne.2) cycle
          write(*,*) 'EDGE CURVE nc = ',nc
        endif
        if (idec_Hermite.eq.1) then
ccc          if (CURVES(nc)%Type.ne.'HermCur') cycle
          if ((CURVES(nc)%Type.ne.'5Bezier').and.
     .        (CURVES(nc)%Type.ne.'7Bezier')) cycle
        endif
        if (idec_rectangle.eq.1) then
          idec=0
          do ir=1,CURVES(nc)%NrFig
            call decode(abs(CURVES(nc)%FigNo(ir)), nr,lab)
            call locate(nr,nrect,nr_rect, ii)
            if (ii.ne.0) idec=1
          enddo
          if (idec.eq.0) cycle
        endif
        if (idec_curve.eq.1) then
          call locate(nc,ncurv,nr_curv, ii)
          if (ii.eq.0) cycle
        endif
        if (idec_bundle.eq.1) then
          idec=0
          do ip=1,2
            if (CURVES(nc)%EndPoNo(ip).eq.npoint) idec=1
          enddo
          if (idec.eq.0) cycle
        endif
c
        nr_curves = nr_curves+1
        do j=1,NRSUB
          ic=ic+1
          COLORN(ic)=ncol
          eta = (j-1)*DX
          call curve(nc,eta, x,dxdeta)
c
c  .......transform into the observer's system
          call trobs(x, xob)
          SEGMENTS(1:2,1,ic) = xob(1:2)
          xmax(1:2) = max(xmax(1:2),xob(1:2))
          xmin(1:2) = min(xmin(1:2),xob(1:2))
          eta = j*DX
          call curve(nc,eta, x,dxdeta)
c
c  .......transform into the observer's system
          call trobs(x, xob)
          SEGMENTS(1:2,2,ic) = xob(1:2)
          xmax(1:2) = max(xmax(1:2),xob(1:2))
          xmin(1:2) = min(xmin(1:2),xob(1:2))
          if (iprint.eq.1) then
            write(*,7003) ic,SEGMENTS(1:2,1,ic),SEGMENTS(1:2,2,ic)
 7003       format('GMPwframe: ic, COORD = ',i4,2(2f8.3,2x))
          endif
        enddo
c
c  .....display the normals
        if (idec_normal.eq.1) then
ccc          if (CURVES(nc)%Type.ne.'HermCur') cycle
          if ((CURVES(nc)%Type.ne.'5Bezier').and.
     .        (CURVES(nc)%Type.ne.'7Bezier')) cycle
          write(*,*) 'displaying normal for nc = ',nc
          do j=1,2
            ic=ic+1
            COLORN(ic)=ncol
            np = CURVES(nc)%EndPoNo(j)
            x(1:3) = POINTS(np)%Rdata(1:3)
            call trobs(x, xob)
            SEGMENTS(1:2,1,ic) = xob(1:2)
            xmax(1:2) = max(xmax(1:2),xob(1:2))
            xmin(1:2) = min(xmin(1:2),xob(1:2))
            x(1:3) = x(1:3) + 2.d0*POINTS(np)%Rdata(4:6)
            call trobs(x, xob)
            SEGMENTS(1:2,2,ic) = xob(1:2)
            xmax(1:2) = max(xmax(1:2),xob(1:2))
            xmin(1:2) = min(xmin(1:2),xob(1:2))
          enddo
        endif
      enddo
      nr_segments=ic
ccc      write(*,*) 'xmin,xmax'
ccc      do i=1,2
ccc        write(*,*) xmin(i),xmax(i)
ccc      enddo
c
c  ...compute object's dimensions and coordinates of its
c     central point
      DIMOB(1) = (xmax(1)-xmin(1))/2.d0
      DIMOB(2) = (xmax(2)-xmin(2))/2.d0
      XCENTR(1) = (xmax(1)+xmin(1))/2.d0
      XCENTR(2) = (xmax(2)+xmin(2))/2.d0
c
      XCIM(1) = XCENTR(1)
      XCIM(2) = XCENTR(2)
c
c  ...decide whether to scale by x or y axes
      q = DIMOB(1)/DIMOB(2)
      q1 = xlength/ylength
      if (q.gt.q1) then
        DIMIM = DIMOB(1)
        SIZE = xlength/2.d0
      else
        DIMIM = DIMOB(2)
        SIZE = ylength/2.d0
      endif
c
c
c***********************************************************************
c
 20   continue
c
c  ...select the window
      call selwin(Iwind)
      write(*,*) '...PLEASE WAIT, PREPARING IMAGE'
c
c  ...display the wireframe
      do ic=1,nr_segments
        ncol = COLORN(ic)
        do k=1,2
          xyp(1:2,k) = SEGMENTS(1:2,k,ic)
          xyp(1:2,k) = (xyp(1:2,k)-XCIM(1:2))/DIMIM*SIZE
     .               + XCWIN(1:2)
        enddo
        call drawline(xyp(1,1),xyp(2,1), xyp(1,2),xyp(2,2), ncol)
      enddo
c
c
c  ...close window
      call dpborder(-1)
      write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
      write(*,*) 'WINDOW TO CONTINUE ...'
      call closwind(iwindnum)
c
      write(*,*) ' SELECT OPTION :'
      write(*,*) ' 0 - EXIT'
      write(*,*) ' 1 - CHANGE THE POINT OF VIEW (WITH RESCALING)'
      write(*,*) ' 2 - CHANGE THE DEGREE OF COARSENESS FOR OVALS'
      write(*,*) ' 3 - DISPLAY ALL CURVES'
      write(*,*) ' 4 - DISPLAY ONLY HERMITE CURVES'
      write(*,*) ' 5 - DISPLAY EDGES OF SELECTED RECTANGLES'
      write(*,*) ' 6 - DISPLAY A POINT BUNDLE'
      write(*,*) ' 7 - DISPLAY CURVES ON SHARP EDGES ONLY'
      write(*,*) ' 8 - DISPLAY NORMALS'
      write(*,*) ' 9 - DISPLAY SELECTED CURVES'
      write(*,*) '10 - CHANGE THE CENTRAL POINT OF THE IMAGE AND',
     .           ' RESCALE THE IMAGE'
      read(*,*) ichoice
c
      select case(ichoice)
      case(0)
        deallocate(SEGMENTS,COLORN)
        return
      case(1)
        write(*,*) 'SET NEW VECTOR NORMAL TO THE PROJECTION PLANE'
        write(*,*) 'psi,theta = ',psi/pi*180,theta/pi*180
        read(*,*)  psi,theta
        psi = psi/180*pi
        theta = theta/180*pi
        RN(1) = cos(psi)*cos(theta)
        RN(2) = cos(psi)*sin(theta)
        RN(3) = sin(psi)
        call cartobs
        xmax(1:2) = bign
        xmin(1:2) = bigp
      case(2)
        write(*,*)'SET NEW DEGREE OF COARSENESS FOR OVALS.'
        write(*,*)'MUST BE BETWEEN 1 AND 5'
        write(*,*)'CURRENT VALUE IS ', nrsub1
        read(*,*) nrsub1
        if (nrsub1.gt.5) then
          write(*,*)'HAVE CHANGED TO 5'
          nrsub1 = 5
        endif
        if (nrsub1.lt.1) then
          write(*,*)'HAVE CHANGED TO 1'
          nrsub1 = 1
        endif
        NRSUB=2**nrsub1
        DX=1.d0/NRSUB
      case(3)
        idec_Hermite=0
        idec_rectangle=0
        idec_curve=0
        idec_bundle=0
        idec_normal=0
        idec_edge=0
      case(4)
c
c  .....reset to the default
        idec_Hermite=1
      case(5)
        idec_rectangle=1
        idec_curve=0
        idec_bundle=0
        idec_edge=0
        write(*,8001) max_nr_rect
 8001   format('GMPwireframe: SET NUMBER OF RECTANGLES',
     .         ' MAX ALLOWED = ',i2)
        read(*,*) nr_rect
        if (nr_rect.gt.max_nr_rect) then
          write(*,8002) max_nr_rect
 8002     format('GMPwireframe: HAVE RESET TO ',i2)
          nr_rect = max_nr_rect
        endif
        write(*,*) 'INPUT THE RECTANGLE NUMBERS'
        read(*,*) nrect(1:nr_rect)
      case(6)
        idec_rectangle=0
        idec_curve=0
        idec_bundle=1
        idec_edge=0
        write(*,*) 'SET THE POINT NUMBER'
        read(*,*) npoint
      case(7)
        idec_rectangle=0
        idec_curve=0
        idec_bundle=0
        idec_edge=1
      case(8)
        idec_normal=1
      case(9)
        idec_rectangle=0
        idec_curve=1
        idec_bundle=0
        idec_edge=0
        write(*,8003) max_nr_curv
 8003   format('GMPwireframe: SET NUMBER OF CURVES',
     .         ' MAX ALLOWED = ',i2)
        read(*,*) nr_curv
        if (nr_curv.gt.max_nr_curv) then
          write(*,8004) max_nr_curv
 8004     format('GMPwireframe: HAVE RESET TO ',i2)
          nr_curv = max_nr_curv
        endif
        write(*,*) 'INPUT THE CURVE NUMBERS'
        read(*,*) ncurv(1:nr_curv)
      case(10)
        write(*,*) 'CLICK AT THE NEW CENTRAL POINT'
        call xmousepos(im(1),im(2))
        XCIM(1) = (im(1) - XCWIN(1))/SIZE*DIMIM + XCIM(1)
        XCIM(2) = (im(2) - XCWIN(2))/SIZE*DIMIM + XCIM(2)
        write(*,*) 'XCIM = ',XCIM
        write(*,*) 'SET MAGNIFICATION FACTOR'
        read(*,*) rnsc
        rnsc = abs(rnsc)
        DIMIM=DIMIM/rnsc
        go to 20
      end select
c
      goto 10
      end subroutine GMPwireframe

#endif
