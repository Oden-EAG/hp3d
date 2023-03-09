!----------------------------------------------------------------------
!> @brief Reconstruction of complete connectivities from partial
!!        connectivities used in the new format
!> @date Mar 2023
!----------------------------------------------------------------------
subroutine connect
!
      use GMP
      use element_data
      implicit none
!
      integer, dimension(:,:), pointer :: iwork
!
!  ...maximum anticipated number of curves attached to a point
      integer, parameter :: maxnrcurv=100
!  ...maximum anticipated number of figures attached to a point
      integer, parameter :: maxnrfig=100
!
      integer :: i1,i2,i3,i4,ic,ii,j,j1,j2,j3,j4,k,k1,k2,lab
      integer :: nc,nc1,nc2,nc3,nc4,nh,nick,nick1,nick2,np,np1,np2
      integer :: nr,nt,npyr,ntet,nrcurv,nrfig,nrfig3,nrfig4,nvoid
!
!  ...an auxiliary list
      integer :: listaux(20)
!
      integer :: iprint=0
      if (iprint.eq.1) then
        write(*,*) 'connect: DEBUGGING...'
      endif
!
!----------------------------------------------------------------------
!
!
!  Step 1: Reconstruct points to curves connectivities
!
      allocate(iwork(1:maxnrcurv,1:NRPOINT), STAT=ic )
      iwork(1:maxnrcurv,1:NRPOINT) = 0
!
!  ...loop through curves
      do nc=1,NRCURVE
!
!  .....loop through the curve endpoints
        do j=1,2
          np = CURVES(nc)%EndPoNo(j)
!
!  .......check if the curve is on the list of connected curves
!         to the point, and add it if it is not
          nrcurv = POINTS(np)%NrCurv
          call locate(nc,iwork(1:nrcurv,np),nrcurv, ii)
          if (ii.eq.0) then
            ii = nrcurv+1
            if (ii.gt.maxnrcurv) then
              write(*,*) 'connect: INCREASE maxnrcurv '
              stop
            endif
            POINTS(np)%NrCurv = ii
            iwork(ii,np) = nc
          endif
        enddo
      enddo
      do np=1,NRPOINT
        allocate(POINTS(np)%CurvNo(1:POINTS(np)%NrCurv))
        do j=1,POINTS(np)%NrCurv
          POINTS(np)%CurvNo(j) = iwork(j,np)
        enddo
        if (iprint.eq.1) then
          write(*,7037) np
 7037     format('connect: CURVES CONNECTED TO POINT ',i6)
          write(*,7038) POINTS(np)%CurvNo(1:POINTS(np)%NrCurv)
 7038     format(10i7)
        endif
      enddo
      if (iprint.eq.1) call pause
      deallocate(iwork)
!
!----------------------------------------------------------------------
!
!  Step 2: Reconstruct rectangle to edge curves connectivities
!
!  ...loop through rectangles
      do nr=1,NRRECTA
!
!  .....loop through rectangle edges
        do j=1,4
!
!  .......local vertex endpoints numbers are
          j1 = QUADR_EDGE_TO_VERT(1,j)
          j2 = QUADR_EDGE_TO_VERT(2,j)
!
!  .......loop through curves connected to the endpoints of the
!         edge and look for a common curve
          np1 = RECTANGLES(nr)%VertNo(j1)
          np2 = RECTANGLES(nr)%VertNo(j2)
          do i1=1,POINTS(np1)%NrCurv
            do i2=1,POINTS(np2)%NrCurv
              if (POINTS(np1)%CurvNo(i1).eq.POINTS(np2)%CurvNo(i2)) then
                RECTANGLES(nr)%EdgeNo(j) = POINTS(np1)%CurvNo(i1)
                goto 20
              endif
            enddo
          enddo
          write(*,7001) nr,j
 7001     format('connect: HAVE NOT FOUND EDGE CURVE FOR RECTANGLE',i5,' AND EDGE ',i1)
          call print_GMP
          stop
 20       continue
!
!  .......determine orientation
          nc = RECTANGLES(nr)%EdgeNo(j)
          if (CURVES(nc)%EndPoNo(1).eq.np2) then
            RECTANGLES(nr)%EdgeNo(j) = - RECTANGLES(nr)%EdgeNo(j)
          elseif (CURVES(nc)%EndPoNo(1).ne.np1) then
            write(*,*) 'connect: INCONSISTENCY 1 '
            write(*,*) 'nc,CURVES(nc)%EndPoNo(1:2) = ',nc,CURVES(nc)%EndPoNo(1:2)
            write(*,*) 'nr,RECTANGLES(nr)%VertNo(1:4) = ',nr,RECTANGLES(nr)%VertNo(1:4)
            call print_GMP
            stop
          endif
!
!  .....end of loop through rectangle edges
        enddo
!
!  ...end of loop through rectangles
      enddo
!
!----------------------------------------------------------------------
!
!  Step 3: Reconstruct triangles to edge curves connectivities
!
!  ...loop through triangles
      do nt=1,NRTRIAN
!
!  .....loop through triangle edges
        do j=1,3
!
!  .......local vertex endpoints numbers are
          j1 = TRIAN_EDGE_TO_VERT(1,j)
          j2 = TRIAN_EDGE_TO_VERT(2,j)
!
!  .......loop through curves connected to the endpoints of the
!         edge and look for a common curve
          np1 = TRIANGLES(nt)%VertNo(j1)
          np2 = TRIANGLES(nt)%VertNo(j2)
          do i1=1,POINTS(np1)%NrCurv
            do i2=1,POINTS(np2)%NrCurv
              if (POINTS(np1)%CurvNo(i1).eq.POINTS(np2)%CurvNo(i2)) then
                TRIANGLES(nt)%EdgeNo(j) = POINTS(np1)%CurvNo(i1)
                goto 50
              endif
            enddo
          enddo
          write(*,7004) nt,j
 7004     format('connect: HAVE NOT FOUND EDGE CURVE FOR TRIANGLE',i8,' AND EDGE ',i1)
          call print_GMP
          stop
 50       continue
!
!  .......determine orientation
          nc = TRIANGLES(nt)%EdgeNo(j)
          if (CURVES(nc)%EndPoNo(1).eq.np2) then
            TRIANGLES(nt)%EdgeNo(j) = - TRIANGLES(nt)%EdgeNo(j)
          elseif (CURVES(nc)%EndPoNo(1).ne.np1) then
            write(*,*) 'connect: INCONSISTENCY 6 '
            stop
          endif
!
!  .....end of loop through triangle edges
        enddo
        if (iprint.eq.1) then
          write(*,7010) nt,TRIANGLES(nt)%EdgeNo(1:3)
 7010     format('connect: nt, EDGES = ',i8,3x,3i8)
        endif
!
!  ...end of loop through triangles
      enddo
!
!----------------------------------------------------------------------
!
!  Step 4: Reconstruct curves to adjacent figures connectivities
!
      allocate(iwork(1:maxnrfig,1:NRCURVE), STAT=ic )
      iwork(1:maxnrfig,1:NRCURVE) = 0
!
!  ...loop through rectangles
      do nr=1,NRRECTA
!
!  .....loop through rectangle edges
        do j=1,4
!
          nc = abs(RECTANGLES(nr)%EdgeNo(j))
!
!  .......check if the rectangle is on the list of figures
!         connected to the curve and add it to the list if it is not
          nrfig = CURVES(nc)%NrFig
          nick = sign(nr*10+2,RECTANGLES(nr)%EdgeNo(j))
          call locate(nick,iwork(1:nrfig,nc),nrfig, ii)
          if (ii.eq.0) then
            nrfig = nrfig+1
            if (nrfig.gt.maxnrfig) then
              write(*,*) 'connect: INCREASE maxnrfig'
              stop
            endif
            CURVES(nc)%NrFig = nrfig
            iwork(nrfig,nc) = nick
          endif
        enddo
      enddo
!
!  ...loop through triangles
      do nt=1,NRTRIAN
!
!  .....loop through triangle edges
        do j=1,3
!
          nc = abs(TRIANGLES(nt)%EdgeNo(j))
!
!  .......check if the triangle is on the list of figures
!         connected to the curve and add it to the list if it is not
          nrfig = CURVES(nc)%NrFig
          nick = sign(nt*10+1,TRIANGLES(nt)%EdgeNo(j))
          call locate(nick,iwork(1:nrfig,nc),nrfig, ii)
          if (ii.eq.0) then
            nrfig = nrfig+1
            if (nrfig.gt.maxnrfig) then
              write(*,*) 'connect: INCREASE maxnrfig'
              stop
            endif
            CURVES(nc)%NrFig = nrfig
            iwork(nrfig,nc) = nick
          endif
        enddo
      enddo
!
      do nc=1,NRCURVE
        allocate(CURVES(nc)%FigNo(1:CURVES(nc)%NrFig), STAT=ic)
        do j=1,CURVES(nc)%NrFig
          CURVES(nc)%FigNo(j) = iwork(j,nc)
        enddo
      enddo
      deallocate(iwork)
      if (NDIM.eq.2) return
!
!----------------------------------------------------------------------
!
!  Step 5: Reconstruct hexa to edge curves and hexa to face
!          rectangles connectivities
!
!  ...loop through hexahedra
      do nh=1,NRHEXAS
!
!  .....loop through the hexa edges
        do j=1,12
!
!  .......the edge endpoints vertex numbers are
          j1 =  BRICK_EDGE_TO_VERT(1,j)
          j2 =  BRICK_EDGE_TO_VERT(2,j)
!
!  .......the corresponding points are
          np1 = HEXAS(nh)%VertNo(j1)
          np2 = HEXAS(nh)%VertNo(j2)
!
!  .......loop through curves connected to the points and look
!         for a common curve
          do i1=1,POINTS(np1)%NrCurv
            do i2=1,POINTS(np2)%NrCurv
              if (POINTS(np1)%CurvNo(i1).eq.POINTS(np2)%CurvNo(i2)) then
                HEXAS(nh)%EdgeNo(j) = POINTS(np1)%CurvNo(i1)
                goto 30
              endif
            enddo
          enddo
          write(*,7002) nh,j
 7002     format('connect: HAVE NOT FOUND EDGE CURVE FOR HEXA ',i5,' AND EDGE ',i2)
          stop
 30       continue
!
!  .......determine orientation
          nc = HEXAS(nh)%EdgeNo(j)
          if (CURVES(nc)%EndPoNo(1).eq.np2) then
            HEXAS(nh)%EdgeNo(j) = - HEXAS(nh)%EdgeNo(j)
          elseif (CURVES(nc)%EndPoNo(1).ne.np1) then
            write(*,*) 'connect: INCONSISTENCY 2'
            stop
          endif
!
!  .....end of loop through edges
        enddo
!
!  .....loop through the hexa faces
        do j=1,6
!
!  .......the face edges are
          j1 = BRICK_FACE_TO_EDGE(1,j)
          j2 = BRICK_FACE_TO_EDGE(2,j)
          j3 = BRICK_FACE_TO_EDGE(3,j)
          j4 = BRICK_FACE_TO_EDGE(4,j)
!
!  .......the corresponding edge curves numbers are
          nc1 = abs(HEXAS(nh)%EdgeNo(j1))
          nc2 = abs(HEXAS(nh)%EdgeNo(j2))
          nc3 = abs(HEXAS(nh)%EdgeNo(j3))
          nc4 = abs(HEXAS(nh)%EdgeNo(j4))
!
!  .......look for a rectangle connected to the first two curves
          do i1=1,CURVES(nc1)%NrFig
            nick1 = abs(CURVES(nc1)%FigNo(i1))
            do i2=1,CURVES(nc2)%NrFig
              nick2 = abs(CURVES(nc2)%FigNo(i2))
              if (nick1.eq.nick2) then
                go to 40
              endif
            enddo
          enddo
          write(*,7003) j,nh
 7003     format('connect: HAVE NOT FOUND ',i1,'TH FACE FOR HEXA',i5)
          stop
 40       continue
!
!  .......check if the rectangle is also connected to the
!         remaining edge curves
          nrfig3 = CURVES(nc3)%NrFig
          listaux(1:nrfig3) = abs(CURVES(nc3)%FigNo(1:nrfig3))
          call locate(nick1,listaux,nrfig3, i3)
          if (i3.eq.0) then
            write(*,*) 'connect: INCONSISTENCY 3'
            stop
          endif
          nrfig4 = CURVES(nc4)%NrFig
          listaux(1:nrfig4) = abs(CURVES(nc4)%FigNo(1:nrfig4))
          call locate(nick1,listaux,nrfig4, i4)
          if (i4.eq.0) then
            write(*,*) 'connect: INCONSISTENCY 4'
            stop
          endif
          call decode(nick1, nr,lab)
!
!  .......determine the orientation of the rectangular face
          k1=0; k2=0
          do k=1,4
            if (abs(RECTANGLES(nr)%EdgeNo(k)).eq.nc1) then
              k1=k
            endif
            if (abs(RECTANGLES(nr)%EdgeNo(k)).eq.nc4) then
              k2=k
            endif
          enddo
          nick = k1*10+k2
          select case(nick)
          case(14)
            HEXAS(nh)%FigNo(j) = nr*10 + 0
          case(43)
            HEXAS(nh)%FigNo(j) = nr*10 + 1
          case(32)
            HEXAS(nh)%FigNo(j) = nr*10 + 2
          case(21)
            HEXAS(nh)%FigNo(j) = nr*10 + 3
          case(41)
            HEXAS(nh)%FigNo(j) = nr*10 + 4
          case(12)
            HEXAS(nh)%FigNo(j) = nr*10 + 5
          case(23)
            HEXAS(nh)%FigNo(j) = nr*10 + 6
          case(34)
            HEXAS(nh)%FigNo(j) = nr*10 + 7
          case default
            write(*,*) 'connect: INCONSISTENCY 5'
            stop
          end select
!
!  .....end of loop through faces
        enddo
!
!  ...end of loop through hexahedra
      enddo
!
!----------------------------------------------------------------------
!
!  Step 6: Reconstruct prisms to edge curves and prisms to face
!          rectangles and triangles connectivities
!
!  ...loop through prisms
      do np=1,NRPRISM
        if (iprint.eq.1) then
          write(*,8001) np
 8001     format('connect: PRISM np = ',i8)
        endif
!
!  .....loop through the prism edges
        do j=1,9
!
!  .......the edge endpoints vertex numbers are
          j1 =  PRISM_EDGE_TO_VERT(1,j)
          j2 =  PRISM_EDGE_TO_VERT(2,j)
!
!  .......the corresponding points are
          np1 = PRISMS(np)%VertNo(j1)
          np2 = PRISMS(np)%VertNo(j2)
!
!  .......loop through curves connected to the points and look
!         for a common curve
          do i1=1,POINTS(np1)%NrCurv
            do i2=1,POINTS(np2)%NrCurv
              if (POINTS(np1)%CurvNo(i1).eq.POINTS(np2)%CurvNo(i2)) then
                PRISMS(np)%EdgeNo(j) = POINTS(np1)%CurvNo(i1)
                goto 60
              endif
            enddo
          enddo
          write(*,7005) np,j
 7005     format('connect: HAVE NOT FOUND EDGE CURVE FOR PRISM ',i5,' AND EDGE ',i2)
          call print_GMP
          stop
 60       continue
!
!  .......determine orientation
          nc = PRISMS(np)%EdgeNo(j)
          if (CURVES(nc)%EndPoNo(1).eq.np2) then
            PRISMS(np)%EdgeNo(j) = - PRISMS(np)%EdgeNo(j)
          elseif (CURVES(nc)%EndPoNo(1).ne.np1) then
            write(*,*) 'connect: INCONSISTENCY 6'
            stop
          endif
!
!  .....end of loop through edges
        enddo
        if (iprint.eq.1) then
          write(*,8002) PRISMS(np)%EdgeNo(1:9)
 8002     format('connect: edges = ',9i8)
        endif
!
!  .....loop through the prism faces
        do j=1,5
!
!  .......the face edges are
          j1 = PRISM_FACE_TO_EDGE(1,j)
          j2 = PRISM_FACE_TO_EDGE(2,j)
          j3 = PRISM_FACE_TO_EDGE(3,j)
          j4 = PRISM_FACE_TO_EDGE(4,j)
!
!  .......the corresponding edge curves numbers are
          nc1 = abs(PRISMS(np)%EdgeNo(j1))
          nc2 = abs(PRISMS(np)%EdgeNo(j2))
          nc3 = abs(PRISMS(np)%EdgeNo(j3))
          nc4 = abs(PRISMS(np)%EdgeNo(j4))
!
!  .......look for a figure connected to the first two curves
          do i1=1,CURVES(nc1)%NrFig
            nick1 = abs(CURVES(nc1)%FigNo(i1))
            do i2=1,CURVES(nc2)%NrFig
              nick2 = abs(CURVES(nc2)%FigNo(i2))
              if (nick1.eq.nick2) then
                go to 70
              endif
            enddo
          enddo
          write(*,7006) j,np
 7006     format('connect: HAVE NOT FOUND ',i1,'TH FACE FOR PRISM',i5)
          stop
 70       continue
!
!  .......check if the figure is also connected to the
!         remaining edge curves
          nrfig3 = CURVES(nc3)%NrFig
          listaux(1:nrfig3) = abs(CURVES(nc3)%FigNo(1:nrfig3))
          call locate(nick1,listaux,nrfig3, i3)
          if (i3.eq.0) then
            write(*,*) 'connect: INCONSISTENCY 7'
            write(*,*) '         np,j = ',np,j
            call print_GMP
            stop
          endif
          if (j.ge.3) then
            nrfig4 = CURVES(nc4)%NrFig
            listaux(1:nrfig4) = abs(CURVES(nc4)%FigNo(1:nrfig4))
            call locate(nick1,listaux,nrfig4, i4)
            if (i4.eq.0) then
              write(*,*) 'connect: INCONSISTENCY 8'
              stop
            endif
          endif
          select case(j)
          case(1,2)
            call decode(nick1, nt,lab)
!
!  .........determine the orientation of the triangular face
            k1=0; k2=0
            do k=1,3
              if (abs(TRIANGLES(nt)%EdgeNo(k)).eq.nc1) then
                k1=k
              endif
              if (abs(TRIANGLES(nt)%EdgeNo(k)).eq.nc2) then
                k2=k
              endif
            enddo
            nick = k1*10+k2
            select case(nick)
            case(12)
            PRISMS(np)%FigNo(j) = nt*10 + 0
            case(31)
              PRISMS(np)%FigNo(j) = nt*10 + 1
            case(23)
              PRISMS(np)%FigNo(j) = nt*10 + 2
            case(32)
              PRISMS(np)%FigNo(j) = nt*10 + 3
            case(13)
              PRISMS(np)%FigNo(j) = nt*10 + 4
            case(21)
              PRISMS(np)%FigNo(j) = nt*10 + 5
            case default
              write(*,*) 'connect: INCONSISTENCY 9'
              stop
            end select
          case(3,4,5)
            call decode(nick1, nr,lab)
!
!  .........determine the orientation of the rectangular face
            k1=0; k2=0
            do k=1,4
              if (abs(RECTANGLES(nr)%EdgeNo(k)).eq.nc1) then
                k1=k
              endif
              if (abs(RECTANGLES(nr)%EdgeNo(k)).eq.nc4) then
                k2=k
              endif
            enddo
            nick = k1*10+k2
            select case(nick)
            case(14)
            PRISMS(np)%FigNo(j) = nr*10 + 0
            case(43)
              PRISMS(np)%FigNo(j) = nr*10 + 1
            case(32)
              PRISMS(np)%FigNo(j) = nr*10 + 2
            case(21)
              PRISMS(np)%FigNo(j) = nr*10 + 3
            case(41)
              PRISMS(np)%FigNo(j) = nr*10 + 4
            case(12)
              PRISMS(np)%FigNo(j) = nr*10 + 5
            case(23)
              PRISMS(np)%FigNo(j) = nr*10 + 6
            case(34)
              PRISMS(np)%FigNo(j) = nr*10 + 7
            case default
              write(*,*) 'connect: INCONSISTENCY 10'
              stop
            end select
          end select
!
!  .....end of loop through faces
        enddo
        if (iprint.eq.1) then
          write(*,8003) PRISMS(np)%FigNo(1:5)
 8003     format('connect: faces = ',9i8)
        endif
!
!  ...end of loop through prisms
      enddo
!
!
!----------------------------------------------------------------------
!
!  Step 7: Reconstruct tetra to edge curves and tetra to face
!          triangles connectivities
!
!  ...loop through tetrahedra
      do ntet=1,NRTETRA
!
!  .....loop through the tet edges
        do j=1,6
!
!  .......the edge endpoints vertex numbers are
          j1 =  TETRA_EDGE_TO_VERT(1,j)
          j2 =  TETRA_EDGE_TO_VERT(2,j)
!
!  .......the corresponding points are
          np1 = TETRAS(ntet)%VertNo(j1)
          np2 = TETRAS(ntet)%VertNo(j2)
!
!  .......loop through curves connected to the points and look
!         for a common curve
          do i1=1,POINTS(np1)%NrCurv
            do i2=1,POINTS(np2)%NrCurv
              if (POINTS(np1)%CurvNo(i1).eq.POINTS(np2)%CurvNo(i2)) then
                TETRAS(ntet)%EdgeNo(j) = POINTS(np1)%CurvNo(i1)
                goto 80
              endif
            enddo
          enddo
          write(*,7020) ntet,j
 7020     format('connect: HAVE NOT FOUND EDGE CURVE FOR TETRA ',i5,' AND EDGE ',i2)
          stop
 80       continue
!
!  .......determine orientation
          nc = TETRAS(ntet)%EdgeNo(j)
          if (CURVES(nc)%EndPoNo(1).eq.np2) then
            TETRAS(ntet)%EdgeNo(j) = - TETRAS(ntet)%EdgeNo(j)
          elseif (CURVES(nc)%EndPoNo(1).ne.np1) then
            write(*,*) 'connect: INCONSISTENCY 11'
            stop
          endif
!
!  .....end of loop through edges
        enddo
!
!  .....loop through the tetra faces
        do j=1,4
!
!  .......the face edges are
          j1 = TETRA_FACE_TO_EDGE(1,j)
          j2 = TETRA_FACE_TO_EDGE(2,j)
          j3 = TETRA_FACE_TO_EDGE(3,j)
!
!  .......the corresponding edge curves numbers are
          nc1 = abs(TETRAS(ntet)%EdgeNo(j1))
          nc2 = abs(TETRAS(ntet)%EdgeNo(j2))
          nc3 = abs(TETRAS(ntet)%EdgeNo(j3))
!
!  .......look for a triangle connected to the first two curves
          do i1=1,CURVES(nc1)%NrFig
            nick1 = abs(CURVES(nc1)%FigNo(i1))
            do i2=1,CURVES(nc2)%NrFig
              nick2 = abs(CURVES(nc2)%FigNo(i2))
              if (nick1.eq.nick2) then
                go to 90
              endif
            enddo
          enddo
          write(*,7023) j,ntet
 7023     format('connect: HAVE NOT FOUND ',i1,'TH FACE FOR TETRA',i5)
          stop
 90       continue
!
!  .......check if the triangle is also connected to the
!         remaining edge curve
          nrfig3 = CURVES(nc3)%NrFig
          listaux(1:nrfig3) = abs(CURVES(nc3)%FigNo(1:nrfig3))
          call locate(nick1,listaux,nrfig3, i3)
          if (i3.eq.0) then
            write(*,*) 'connect: INCONSISTENCY 12'
            write(*,*) 'ntet,j = ',ntet,j
            call print_GMP
            stop
          endif
          call decode(nick1, nt,lab)
!
!  .......determine the orientation of the triangular face
          k1=0; k2=0
          do k=1,3
            if (abs(TRIANGLES(nt)%EdgeNo(k)).eq.nc1) then
              k1=k
            endif
            if (abs(TRIANGLES(nt)%EdgeNo(k)).eq.nc2) then
              k2=k
            endif
          enddo
          nick = k1*10+k2
          select case(nick)
          case(12)
            TETRAS(ntet)%FigNo(j) = nt*10 + 0
          case(31)
            TETRAS(ntet)%FigNo(j) = nt*10 + 1
          case(23)
            TETRAS(ntet)%FigNo(j) = nt*10 + 2
          case(32)
            TETRAS(ntet)%FigNo(j) = nt*10 + 3
          case(13)
            TETRAS(ntet)%FigNo(j) = nt*10 + 4
          case(21)
            TETRAS(ntet)%FigNo(j) = nt*10 + 5
          case default
            write(*,*) 'connect: INCONSISTENCY 5'
            stop
          end select
!
!  .....end of loop through faces
        enddo
!
!  ...end of loop through tetrahedra
      enddo
!
!
!----------------------------------------------------------------------
!
!  Step 8: Reconstruct pyramid to edge curves and pyramid to face
!          rectangle and triangles connectivities
!
!  ...loop through pyramids
      do npyr=1,NRPYRAM
        if (iprint.eq.1) then
          write(*,7101) npyr
 7101     format('connect: PYRAMID npyr = ',i8)
        endif
!
!  .....loop through the pyramid edges
        do j=1,8
!
!  .......the edge endpoints vertex numbers are
          j1 =  PYRAM_EDGE_TO_VERT(1,j)
          j2 =  PYRAM_EDGE_TO_VERT(2,j)
!
!  .......the corresponding points are
          np1 = PYRAMIDS(npyr)%VertNo(j1)
          np2 = PYRAMIDS(npyr)%VertNo(j2)
!
!  .......loop through curves connected to the points and look
!         for a common curve
          do i1=1,POINTS(np1)%NrCurv
            do i2=1,POINTS(np2)%NrCurv
              if (POINTS(np1)%CurvNo(i1).eq.POINTS(np2)%CurvNo(i2)) then
                PYRAMIDS(npyr)%EdgeNo(j) = POINTS(np1)%CurvNo(i1)
                goto 160
              endif
            enddo
          enddo
          write(*,7105) npyr,j
 7105     format('connect: HAVE NOT FOUND EDGE CURVE FOR PYRAMID ',i5,' AND EDGE ',i2)
          stop
 160      continue
!
!  .......determine orientation
          nc = PYRAMIDS(npyr)%EdgeNo(j)
          if (CURVES(nc)%EndPoNo(1).eq.np2) then
            PYRAMIDS(npyr)%EdgeNo(j) = - PYRAMIDS(npyr)%EdgeNo(j)
          elseif (CURVES(nc)%EndPoNo(1).ne.np1) then
            write(*,*) 'connect: INCONSISTENCY 16'
            write(*,*) 'np1,np2,nc = ',np1,np2,nc
            write(*,*) 'npyr,j = ',npyr,j
            call print_GMP
            stop
          endif
!
!  .....end of loop through edges
        enddo
        if (iprint.eq.1) then
          write(*,8002) PYRAMIDS(npyr)%EdgeNo(1:8)
        endif
!
!  .....loop through the pyramid faces
        do j=1,5
!
!  .......the face edges are
          j1 = PYRAM_FACE_TO_EDGE(1,j)
          j2 = PYRAM_FACE_TO_EDGE(2,j)
          j3 = PYRAM_FACE_TO_EDGE(3,j)
          j4 = PYRAM_FACE_TO_EDGE(4,j)
!
!  .......the corresponding edge curves numbers are
          nc1 = abs(PYRAMIDS(npyr)%EdgeNo(j1))
          nc2 = abs(PYRAMIDS(npyr)%EdgeNo(j2))
          nc3 = abs(PYRAMIDS(npyr)%EdgeNo(j3))
          nc4 = abs(PYRAMIDS(npyr)%EdgeNo(j4))
!
!  .......look for a figure connected to the first two curves
          do i1=1,CURVES(nc1)%NrFig
            nick1 = abs(CURVES(nc1)%FigNo(i1))
            do i2=1,CURVES(nc2)%NrFig
              nick2 = abs(CURVES(nc2)%FigNo(i2))
              if (nick1.eq.nick2) then
                go to 170
              endif
            enddo
          enddo
          write(*,7106) j,npyr
 7106     format('connect: HAVE NOT FOUND ',i1,'TH FACE FOR PYRAMID',i5)
          stop
 170      continue
!
!  .......check if the figure is also connected to the
!         remaining edge curves
          nrfig3 = CURVES(nc3)%NrFig
          listaux(1:nrfig3) = abs(CURVES(nc3)%FigNo(1:nrfig3))
          call locate(nick1,listaux,nrfig3, i3)
          if (i3.eq.0) then
            write(*,*) 'connect: INCONSISTENCY 17'
            write(*,*) '         np,j = ',np,j
            call print_GMP
            stop
          endif
          if (j.eq.1) then
            nrfig4 = CURVES(nc4)%NrFig
            listaux(1:nrfig4) = abs(CURVES(nc4)%FigNo(1:nrfig4))
            call locate(nick1,listaux,nrfig4, i4)
            if (i4.eq.0) then
              write(*,*) 'connect: INCONSISTENCY 18'
              stop
            endif
          endif
!
          select case(j)
!
!  .......determine the orientation of the rectangular face
          case(1)
            call decode(nick1, nr,lab)
!
!  .........determine the orientation of the rectangular face
            k1=0; k2=0
            do k=1,4
              if (abs(RECTANGLES(nr)%EdgeNo(k)).eq.nc1) then
                k1=k
              endif
              if (abs(RECTANGLES(nr)%EdgeNo(k)).eq.nc4) then
                k2=k
              endif
            enddo
            nick = k1*10+k2
            select case(nick)
            case(14)
            PYRAMIDS(npyr)%FigNo(j) = nr*10 + 0
            case(43)
              PYRAMIDS(npyr)%FigNo(j) = nr*10 + 1
            case(32)
              PYRAMIDS(npyr)%FigNo(j) = nr*10 + 2
            case(21)
              PYRAMIDS(npyr)%FigNo(j) = nr*10 + 3
            case(41)
              PYRAMIDS(npyr)%FigNo(j) = nr*10 + 4
            case(12)
              PYRAMIDS(npyr)%FigNo(j) = nr*10 + 5
            case(23)
              PYRAMIDS(npyr)%FigNo(j) = nr*10 + 6
            case(34)
              PYRAMIDS(npyr)%FigNo(j) = nr*10 + 7
            case default
              write(*,*) 'connect: INCONSISTENCY 110'
              stop
            end select
          case(2,3,4,5)
            call decode(nick1, nt,lab)
!
!  .........determine the orientation of the triangular face
            k1=0; k2=0
            do k=1,3
              if (abs(TRIANGLES(nt)%EdgeNo(k)).eq.nc1) then
                k1=k
              endif
              if (abs(TRIANGLES(nt)%EdgeNo(k)).eq.nc2) then
                k2=k
              endif
            enddo
            nick = k1*10+k2
            select case(nick)
            case(12)
            PYRAMIDS(npyr)%FigNo(j) = nt*10 + 0
            case(31)
              PYRAMIDS(npyr)%FigNo(j) = nt*10 + 1
            case(23)
              PYRAMIDS(npyr)%FigNo(j) = nt*10 + 2
            case(32)
              PYRAMIDS(npyr)%FigNo(j) = nt*10 + 3
            case(13)
              PYRAMIDS(npyr)%FigNo(j) = nt*10 + 4
            case(21)
              PYRAMIDS(npyr)%FigNo(j) = nt*10 + 5
            case default
              write(*,*) 'connect: INCONSISTENCY 19'
              stop
            end select
          end select
!
!  .....end of loop through faces
        enddo
        if (iprint.eq.1) then
          write(*,8003) PYRAMIDS(npyr)%FigNo(1:5)
        endif
!
!  ...end of loop through pyramids
      enddo
!
!----------------------------------------------------------------------
!
!  Step 8: Reconstruct rectangles to blocks connectivities
!
!  ...loop through hexa
      do nh=1,NRHEXAS
        nick = nh*10+2
!
!  .....loop through hexa faces
        do j=1,6
!
          call decode(HEXAS(nh)%FigNo(j), nr,nvoid)
!
!  .......check if the hexa is on the list of blocks
!         connected to the recta and add it to the list if it is not
          call locate(nick,RECTANGLES(nr)%BlockNo,2, ii)
          if (ii.eq.0) then
            if (RECTANGLES(nr)%BlockNo(1).eq.0) then
              RECTANGLES(nr)%BlockNo(1) = nick
            else
              RECTANGLES(nr)%BlockNo(2) = nick
            endif
          endif
        enddo
!
!  ...end of loop through hexa
      enddo
      if (iprint.eq.1) then
        write(*,*) 'connect: HAVE CONNECTED HEXAS TO RECTANGLES'
      endif
!
!  ...loop through prisms
      do np=1,NRPRISM
        nick = np*10+1
!
!  .....loop through prism faces
        do j=3,5
!
          call decode(PRISMS(np)%FigNo(j), nr,nvoid)
!
!  .......check if the prism is on the list of blocks
!         connected to the recta and add it to the list if it is not
          call locate(nick,RECTANGLES(nr)%BlockNo,2, ii)
          if (ii.eq.0) then
            if (RECTANGLES(nr)%BlockNo(1).eq.0) then
              RECTANGLES(nr)%BlockNo(1) = nick
            else
              RECTANGLES(nr)%BlockNo(2) = nick
            endif
          endif
        enddo
!
!  ...end of loop through prisms
      enddo
      if (iprint.eq.1) then
        write(*,*) 'connect: HAVE CONNECTED PRISMS TO RECTANGLES'
      endif
!
!  ...loop through pyramids
      do npyr=1,NRPYRAM
        nick = npyr*10+4
!
!  .....bottom face only
        call decode(PYRAMIDS(npyr)%FigNo(1), nr,nvoid)
!
!  .....check if the pyramid is on the list of blocks
!       connected to the recta and add it to the list if it is not
        call locate(nick,RECTANGLES(nr)%BlockNo,2, ii)
        if (ii.eq.0) then
          if (RECTANGLES(nr)%BlockNo(1).eq.0) then
            RECTANGLES(nr)%BlockNo(1) = nick
          else
            RECTANGLES(nr)%BlockNo(2) = nick
          endif
        endif
!
!  ...end of loop through pyramids
      enddo
      if (iprint.eq.1) then
        write(*,*) 'connect: HAVE CONNECTED PYRAMIDS TO RECTANGLES'
      endif
!
!
!----------------------------------------------------------------------
!
!  Step 9: Reconstruct triangles to blocks connectivities
!
!
!  ...loop through prisms
      do np=1,NRPRISM
        nick = np*10+1
!
!  .....loop through prism bases
        do j=1,2
!
          call decode(PRISMS(np)%FigNo(j), nt,nvoid)
!
!  .......check if the prism is on the list of blocks
!         connected to the triangle and add it to the list if it is not
          call locate(nick,TRIANGLES(nt)%BlockNo,2, ii)
          if (ii.eq.0) then
            if (TRIANGLES(nt)%BlockNo(1).eq.0) then
              TRIANGLES(nt)%BlockNo(1) = nick
            else
              TRIANGLES(nt)%BlockNo(2) = nick
            endif
          endif
        enddo
!
!  ...end of loop through prisms
      enddo
      if (iprint.eq.1) then
        write(*,*) 'connect: HAVE CONNECTED PRISMS TO TRIANGLES'
      endif
!
!  ...loop through tets
      do ntet=1,NRTETRA
        nick = ntet*10+3
!
!  .....loop through the tetra faces
        do j=1,4
!
          call decode(TETRAS(ntet)%FigNo(j), nt,nvoid)
!
!  .......check if the tetrahedron is on the list of blocks
!         connected to the triangle and add it to the list if it is not
          call locate(nick,TRIANGLES(nt)%BlockNo,2, ii)
          if (ii.eq.0) then
            if (TRIANGLES(nt)%BlockNo(1).eq.0) then
              TRIANGLES(nt)%BlockNo(1) = nick
            else
              TRIANGLES(nt)%BlockNo(2) = nick
            endif
          endif
        enddo
!
!  ...end of loop through tetras
      enddo
      if (iprint.eq.1) then
        write(*,*) 'connect: HAVE CONNECTED TETS TO TRIANGLES'
      endif
!
!  ...loop through pyramids
      do npyr=1,NRPYRAM
        nick = npyr*10+4
!
!  .....loop through the pyramid lateral faces
        do j=2,5
!
          call decode(PYRAMIDS(npyr)%FigNo(j), nt,nvoid)
!
!  .......check if the pyramid is on the list of blocks
!         connected to the triangle and add it to the list if it is not
          call locate(nick,TRIANGLES(nt)%BlockNo,2, ii)
          if (ii.eq.0) then
            if (TRIANGLES(nt)%BlockNo(1).eq.0) then
              TRIANGLES(nt)%BlockNo(1) = nick
            else
              TRIANGLES(nt)%BlockNo(2) = nick
            endif
          endif
        enddo
!
!  ...end of loop through pyramids
      enddo
      if (iprint.eq.1) then
        write(*,*) 'connect: HAVE CONNECTED PYRAMIDS TO TRIANGLES'
      endif
!
!
end subroutine connect
!
!
!----------------------------------------------------------------------
!> @brief Deletes connectivity info for GMP objects
!> @date Mar 2023
!----------------------------------------------------------------------
subroutine clean_GMP
!
      use GMP
      implicit none
!
      integer :: np,nc,nt,nr,nh,ntet,npri,npyr
!
      do np=1,NRPOINT
        POINTS(np)%NrCurv=0
        deallocate(POINTS(np)%CurvNo)
      enddo
!
      do nc=1,NRCURVE
        CURVES(nc)%Nrfig=0
        deallocate(CURVES(nc)%FigNo)
      enddo
!
      do nt=1,NRTRIAN
        TRIANGLES(nt)%EdgeNo=0
        TRIANGLES(nt)%BlockNo=0
      enddo
!
      do nr=1,NRRECTA
        RECTANGLES(nr)%EdgeNo=0
        RECTANGLES(nr)%BlockNo=0
      enddo
!
      do nh=1,NRHEXAS
        HEXAS(nh)%EdgeNo=0
        HEXAS(nh)%FigNo=0
      enddo
!
      do ntet=1,NRTETRA
        TETRAS(ntet)%EdgeNo=0
        TETRAS(ntet)%FigNo=0
      enddo
!
      do npri=1,NRPRISM
        PRISMS(npri)%EdgeNo=0
        PRISMS(npri)%FigNo=0
      enddo
!
      do npyr=1,NRPYRAM
        PYRAMIDS(npyr)%EdgeNo=0
        PYRAMIDS(npyr)%FigNo=0
      enddo
!
!
end subroutine clean_GMP
