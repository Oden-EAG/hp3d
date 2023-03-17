c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : this is a LEGACY routine. Do not use it, unless
c                    you revise it!
c//////////////////////////////////////////////////////////////////////
c
c
c----------------------------------------------------------------------
c
c   routine name       - input_GMPdata
c
c----------------------------------------------------------------------
c
c   latest revision    - Aug 09
c
c   purpose            - routine reads in geometry data for a hexa only
c                        GMP model in the COMPACT FORMAT, allocates
c                        dynamically necessary data structure
c                        arrays, and produces automatically all
c                        necessary info for the GMP data structure
c
c
c   arguments in       - none
c
c---------------------------------------------------------------------
c
      subroutine input_GMPdata
c
      use GMP
      use element_data
#include "syscom.blk"
c
      integer, dimension(:,:), pointer:: point_to_blocks
      integer, dimension(:),   pointer:: point_nrbl
      integer, dimension(:,:), pointer:: hexa_to_surf
c
      parameter (mphex = 10)
c
      dimension list(2,20)
c     for i-th block on the list
c       list(1,i) = block connected to an edge
c       list(2,i) = the corresponding edge number
c     or
c       list(1,i) = block connected to a face
c       list(2,i) = the corresponding face number
c
c---------------------------------------------------------------------
c
      iprint=1
c
c  ...read in the dimension of the problem and the manifold
      read(KIN,*)  NDIM,MANDIM
c
c  ...allocate memory for GMP data structure
      call alloc_GMP
c
c---------------------------------------------------------------------
c
c  Step 1: input surfaces
c
      call input_SURFACES(KIN)
c
c----------------------------------------------------------------------
c
c  Step 2: Input points
c
      call input_POINTS(KIN)
c
c----------------------------------------------------------------------
c
c  Step 3: Input hexas
c
c  ...read in number of hexas
      read(KIN,*) NRHEXAS
c
      allocate(hexa_to_surf(6,NRHEXAS))
c
      if (iprint.eq.1) then
        write(*,*) 'input_GMPdata: READING HEXAS...'
      endif
c
      do ih=1,NRHEXAS
        HEXAS(ih)%Type = 'TriLiHex'
        read(KIN,*) HEXAS(ih)%VertNo(1:8)
        read(KIN,*) hexa_to_surf(1:6,ih)
c
c  ...end of loop through hexahedra
      enddo
c
c  ...no other blocks but bricks for a while
      NRPRISM = 0; NRTETRA = 0; NRPYRAM = 0
c
c----------------------------------------------------------------------
c
c  Step 4: Determine point to hexas connectivities
c
      allocate(point_to_blocks(1:mphex,1:NRPOINT))
      allocate(point_nrbl(1:NRPOINT))
      point_nrbl(1:NRPOINT) = 0
c
c  ...loop through hexa
      do nh=1,NRHEXAS
c
        nick = nh*10+2
c
c  .....loop through hexa vertex points
        do jv=1,8
c
          np = HEXAS(nh)%VertNo(jv)
          nrbl = point_nrbl(np)
c
c  .......check if the hexa is on the list of blocks
c         connected to the point and add it to the list if it is not
          call locate(nick,point_to_blocks(1:nrbl,np),nrbl, ii)
          if (ii.eq.0) then
            point_to_blocks(nrbl+1,np) = nick
            point_nrbl(np) = nrbl + 1
          endif
c
c  .....end loop through vertices
        enddo
c
c  ...end of loop through hexa
      enddo
      if (iprint.eq.1) then
        do np=1,NRPOINT
          write(*,7005) np, (point_to_blocks(ii,np),ii=1,point_nrbl(np))
 7005     format('input_GMPdata: np = ',i4,' CONNECTED BLOCKS = ',10i4)
        enddo
        call pause
      endif
c
c----------------------------------------------------------------------
c
c  Step 5a: Count curves and determine temporary connectivities
c
      NRCURVE=0
c
      do nh=1,NRHEXAS
        HEXAS(nh)%EdgeNo(1:12) = 0
      enddo
c
c  ...loop through hexa
      do nh=1,NRHEXAS
c
        nick = nh*10+2
c
c  .....loop through the hexa edge curves
        do je=1,12
c
c  .......if the edge has not been connected yet
          if (HEXAS(nh)%EdgeNo(je).eq.0) then
c
c  .........initiate the list of blocks connected to the edge
c           with the current block
            ibl=0
c
c  .........determine the endpoints of the edge
            iv1 = BRICK_EDGE_TO_VERT(1,je)
            np1 = HEXAS(nh)%VertNo(iv1)
            iv2 = BRICK_EDGE_TO_VERT(2,je)
            np2 = HEXAS(nh)%VertNo(iv2)
c
c  .........loop through the hexas connected to the first endpoint
            nrbl = point_nrbl(np1)
            do k=1,nrbl
              nick1 = point_to_blocks(k,np1)
              call decode(nick1, nh1,lab)
              select case(lab)
c
c  ...........hexahedron
              case(2)
c
c  .............locate the vertex coinciding with the point
                call locate(np1,HEXAS(nh1)%VertNo(1:8),8, kv1)
c
c  .............loop through the edges sharing the vertex
                iflag=0
                do l=1,3
                  ie = BRICK_VERT_TO_EDGE(l,kv1)
c
c  ...............look for the other vertex point
                  do ive=1,2
                    iv3 = BRICK_EDGE_TO_VERT(ive,ie)
                    np3 = HEXAS(nh1)%VertNo(iv3)
                    if (np2.eq.np3) iflag=1
                  enddo
                  if (iflag.eq.1) exit
c
c  .............end of loop through the edge sharing the vertex
                enddo
                if (iflag.eq.0) go to 10
                nc = HEXAS(nh1)%EdgeNo(ie)
c
c  .............if the curve has already been counted
                if (nc.ne.0) then
c
c  ...............store the connection for the original hexa
                  HEXAS(nh)%EdgeNo(je) = nc
                  go to 20
                else
c
c  ...............add the block to the list of blocks
c                 to be connected to the curve
                  ibl=ibl+1
                  list(1,ibl) = nick1; list(2,ibl) = ie
c
c  ...............proceed to the next connected block
                  go to 10
                endif
              case default
                write(*,*) 'input_GMPdata: WRONG lab = ',lab
                stop 1
              end select
c
 10           continue
c
c  .........end of loop through connected blocks
            enddo
c
c  .........no block connected to the edge has been found, the
c           edge curve has to be created
            NRCURVE = NRCURVE+1
            nc = NRCURVE
            if (iprint.eq.1) then
              write(*,7003) nc,nh,je
 7003         format('input_GMPdata: HAVE DEFINED CURVE ',i4,
     .               ' FOR HEXA ', i4,' AND EDGE ',i2)
              write(*,*) 'CONNECTED BLOCKS WITH EDGE NUMBERS'
              write(*,7007) (list(1,jj)/10,jj=1,ibl)
 7007         format(10i5)
              write(*,7007) (list(2,jj),jj=1,ibl)
              call pause
            endif
c
c  .........connect all blocks from the list to the curve
            nrbl = ibl
            do ibl=1,nrbl
              call decode(list(1,ibl), nn,lab)
              select case(lab)
              case(2)
                ie = list(2,ibl)
                HEXAS(nn)%EdgeNo(ie) = nc
              end select
            enddo
c
c  .......if the edge has not been connected yet
          endif
 20       continue
c
c  .....end of loop through edges of the hexa
        enddo
c
c  ...end of loop through hexas
      enddo
c
c  ...define curves..................
      do nc=1,NRCURVE
        CURVES(nc)%Type = 'void'
      enddo
c
c  ...loop through hexa
      do nh=1,NRHEXAS
c
c  .....loop through the hexa edge curves
        do je=1,12
c
c  .......get the connected curve
          nc = HEXAS(nh)%EdgeNo(je)
c
c  .......if the curve is visited for the first type
          if (CURVES(nc)%Type.eq.'void') then
c
c  .........determine the endpoints of the edge
            iv1 = BRICK_EDGE_TO_VERT(1,je)
            np1 = HEXAS(nh)%VertNo(iv1)
            iv2 = BRICK_EDGE_TO_VERT(2,je)
            np2 = HEXAS(nh)%VertNo(iv2)
            CURVES(nc)%Type = 'Seglin'
            CURVES(nc)%EndPoNo(1) = np1
            CURVES(nc)%EndPoNo(2) = np2
            CURVES(nc)%NrFig = 0
          endif
c
c  .....end of loop through edges of the hexa
        enddo
c
c  ...end of loop through hexas
      enddo
c
c  ...erase the temporary block to edge curves connectivities
      do nh=1,NRHEXAS
        HEXAS(nh)%EdgeNo(1:12) = 0
      enddo
      if (iprint.eq.1) then
        do nc = 1,NRCURVE
          write(*,7001) nc,(CURVES(nc)%EndPoNo(ii),ii=1,2)
 7001     format('input_GMPdata: nc,points = ',i4,2x,2i5)
        enddo
        call pause
      endif
c
c----------------------------------------------------------------------
c
c  ...no triangles for a while
      NRTRIAN = 0
c
c  Step 6: count rectas and determine temporary hexas to face rectas
c          connectivities
c
      NRRECTA=0
c
      do nh=1,NRHEXAS
        HEXAS(nh)%FigNo(1:6) = 0
      enddo
c
c  ...loop through hexa
      do nh=1,NRHEXAS
c
        nick = nh*10+2
c
c  .....loop through the hexa face rectangles
        do jf=1,6
c
c  .......if the face has not been connected yet
          if (HEXAS(nh)%FigNo(jf).eq.0) then
c
c  .........initiate the list of blocks connected to the face
c           withe the current block
            ibl=0
c
c  .........determine vertex points on the face
            iv = BRICK_FACE_TO_VERT(1,jf)
            np1 = HEXAS(nh)%VertNo(iv)
            iv = BRICK_FACE_TO_VERT(2,jf)
            np2 = HEXAS(nh)%VertNo(iv)
            iv = BRICK_FACE_TO_VERT(3,jf)
            np3 = HEXAS(nh)%VertNo(iv)
            iv = BRICK_FACE_TO_VERT(4,jf)
            np4 = HEXAS(nh)%VertNo(iv)
c
c  .........loop through the hexas connected to the point
            nrbl = point_nrbl(np1)
            do k=1,nrbl
              nick1 = point_to_blocks(k,np1)
              call decode(nick1, nh1,lab)
              select case(lab)
c
c  ...........hexahedron
              case(2)
c
c  .............locate the vertex coinciding with the point
                call locate(np1,HEXAS(nh1)%VertNo(1:8),8, iv1)
c
c  .............loop through the faces sharing the point
                do l=1,3
                  if = BRICK_VERT_TO_FACE(l,iv1)
c
c  ...............look for two two other common points
                  iflag=0
                  do ivf=1,4
                    iv = BRICK_FACE_TO_VERT(ivf,if)
                    np = HEXAS(nh1)%VertNo(iv)
                    if (np.eq.np2) iflag=iflag+1
                    if (np.eq.np3) iflag=iflag+1
                  enddo
                  if (iflag.eq.2) exit
                enddo
                if (iflag.ne.2) go to 30
c
c  .............look for a connected rectangle
                nf = HEXAS(nh1)%FigNo(if)
c
c  .............if the rectangle has already been created
                if (nf.ne.0) then
c
c  ...............store the connection for the original hexa
                  HEXAS(nh)%FigNo(jf) = nf
                  go to 40
                else
c
c  ...............add the block to the list of blocks
c                 to be connected to the curve
                  ibl=ibl+1
                  list(1,ibl) = nick1; list(2,ibl) = if
c
c  ...............proceed to the next connected block
                  go to 30
                endif
              case default
                write(*,*) 'input_GMPdata: WRONG lab = ',lab
                stop 1
              end select
c
 30           continue
c
c  .........end of loop through connected blocks
            enddo
c
c  .........no block connected to the face has been found, the
c           face rectangle has to be created
            NRRECTA = NRRECTA + 1
            nr = NRRECTA
            if (iprint.eq.1) then
              write(*,7008) nr,nh,jf
 7008         format('input_GMPdata: HAVE DEFINED RECTA ',i4,
     .               ' FOR HEXA ', i4,' AND FACE ',i2)
              write(*,*) 'CONNECTED BLOCKS WITH FACE NUMBERS'
              write(*,7007) (list(1,jj)/10,jj=1,ibl)
              write(*,7007) (list(2,jj),jj=1,ibl)
              call pause
            endif
c
c  .........connect all blocks from the list to the rectangle
            nrbl = ibl
            do ibl=1,nrbl
              call decode(list(1,ibl), nn,lab)
              select case(lab)
              case(2)
                if = list(2,ibl)
                HEXAS(nn)%FigNo(if) = nr
              end select
            enddo
c
c  .......if the face has not been connected yet
          endif
 40       continue
c
c  .....end of loop through faces of the hexa
        enddo
c
c  ...end of loop through hexas
      enddo
c
c
c  ...define the rectangles......................
      do nr=1,NRRECTA
        RECTANGLES(nr)%Type = 'void'
      enddo
c
c  ...loop through hexa
      do nh=1,NRHEXAS
c
c  .....loop through the hexa face rectangles
        do jf=1,6
          nr = HEXAS(nh)%FigNo(jf)
c
c  .......if visited for the first time
          if (RECTANGLES(nr)%Type.eq.'void') then
c
c  .........determine vertex points on the face
            iv = BRICK_FACE_TO_VERT(1,jf)
            np1 = HEXAS(nh)%VertNo(iv)
            iv = BRICK_FACE_TO_VERT(2,jf)
            np2 = HEXAS(nh)%VertNo(iv)
            iv = BRICK_FACE_TO_VERT(3,jf)
            np3 = HEXAS(nh)%VertNo(iv)
            iv = BRICK_FACE_TO_VERT(4,jf)
            np4 = HEXAS(nh)%VertNo(iv)
            RECTANGLES(nr)%Type = 'BilQua'
            RECTANGLES(nr)%VertNo(1) = np1
            RECTANGLES(nr)%VertNo(2) = np2
            RECTANGLES(nr)%VertNo(3) = np3
            RECTANGLES(nr)%VertNo(4) = np4
            RECTANGLES(nr)%BlockNo(1:2) = 0
c
          endif
c
c  .....end of loop through faces of the hexa
        enddo
c
c  ...end of loop through hexas
      enddo
c
c  ...erase the temporary block to face rectangles connectivities
      do nh=1,NRHEXAS
        HEXAS(nh)%FigNo(1:6) = 0
      enddo
      if (iprint.eq.1) then
        do nr = 1,NRRECTA
          write(*,7002) nr,(RECTANGLES(nr)%VertNo(ii),ii=1,4)
 7002     format('input_GMPdata: nr,points = ',i4,2x,4i5)
        enddo
        call pause
      endif
c
c----------------------------------------------------------------------
c
c  Step 7: complete connectivities
c
      deallocate(point_to_blocks,point_nrbl)
      call connect
c
c
c----------------------------------------------------------------------
c
c  Step 8: define implicit curves and rectangles
c
      call def_impl(hexa_to_surf)
c
c
c  ...clean up...
      deallocate (hexa_to_surf)
c
c
      end subroutine input_GMPdata
c
c
c
!----------------------------------------------------------------------
!   latest revision    - Nov 08
!
!   purpose            - Routine inputs data on surfaces
!
!   arguments
!     in:
!              Input   - input file
!---------------------------------------------------------------------
!
      subroutine input_SURFACES(Input)
!
      use control
      use GMP
#include "syscom.blk"
!
!---------------------------------------------------------------------
!
      iprint=0
!
!  ...read in number of surfaces
      read(Input,*) NRSURFS
      write(*,*)'input_DEFAULT: READING SURFACES...'
      write(*,*)'           NRSURFS = ',NRSURFS
!
      if (MAXSU.lt.NRSURFS) then
        write(*,7001)
 7001   format('input_SURFACES: INCREASE MAXSU')
        stop
      endif
      if (iprint.eq.1) then
        write(*,7002) NRSURFS
 7002   format(' input_SURFACES: NRSURFS = ',i5,' READING SURFACES...')
      endif
!
      do ns=1,NRSURFS
        read(Input,*) SURFACES(ns)%Type
        if (iprint.eq.1) then
          write(*,7003) ns,SURFACES(ns)%Type
 7003     format('input_SURFACES: ns = ',i3,' SURFACE TYPE = ',a10)
        endif
!
        select case(SURFACES(ns)%Type)
!
!  .......a plane normal to a given vector and passing through
!         a point, orientation specified by the normal
          case('VecPt')
            allocate(SURFACES(ns)%Rdata(6))
            read(Input,*) SURFACES(ns)%Rdata(1:3)
            read(Input,*) SURFACES(ns)%Rdata(4:6)
!
!  .......plane passing through three given points, A,B,C
!         orientation specified by AB x AC
          case('ThrPt')
            allocate(SURFACES(ns)%Rdata(9))
            read(Input,*) SURFACES(ns)%Rdata(1:3)
            read(Input,*) SURFACES(ns)%Rdata(4:6)
            read(Input,*) SURFACES(ns)%Rdata(7:9)
!
!  .......a sphere with radius and central point
          case('Sphere')
            allocate(SURFACES(ns)%Rdata(4))
            read(Input,*) SURFACES(ns)%Rdata(1:3)
            read(Input,*) SURFACES(ns)%Rdata(4)
!
!    .....a cylinder defined with an endpoint of its axis, direction
!         vector and radius
          case('Cylinder')
            allocate(SURFACES(ns)%Rdata(7))
            read(Input,*) SURFACES(ns)%Rdata(1:3)
            read(Input,*) SURFACES(ns)%Rdata(4:6)
            read(Input,*) SURFACES(ns)%Rdata(7)
!
!  .......ellipsoid with semi-axes a,b,c along x,y,z and given central
!         point
          case('Ellipsoid')
            allocate(SURFACES(ns)%Rdata(6))
!
!  .........read coordinates of the central point
            read(Input,*) SURFACES(ns)%Rdata(1:3)
!
!  .........read semiaxes lengths
            read(Input,*) SURFACES(ns)%Rdata(4:6)
!
!  .......cone
          case('Cone')
            allocate(SURFACES(ns)%Rdata(7))
!
!  .........read coordinates of the central point
            read(Input,*) SURFACES(ns)%Rdata(1:3)
!
!  .........read direction vector
            read(Input,*) SURFACES(ns)%Rdata(4:6)
!
!  .........read the cone angle
            read(Input,*) SURFACES(ns)%Rdata(7)
!
!  .......reconstructed surface
          case('RecSurf')
!
!  .......other cases
          case default
            write(*,7004)
 7004       format('input_SURFACES: UNKNOWN SURFACE TYPE')
            write(*,*) SURFACES(ns)%Type
            stop
        endselect
      enddo
!
!
      end subroutine input_SURFACES
!
!
!
!----------------------------------------------------------------------
!   latest revision    - Nov 10
!
!   purpose            - Routine inputs data on GMP points
!
!   arguments
!     in:
!              Input   - input file
!---------------------------------------------------------------------
!
      subroutine input_POINTS(Input)
!
      use control
      use GMP
#include "syscom.blk"
!
      character(10) type
!
!---------------------------------------------------------------------
!
      iprint=0
!
!
!  ...read in number of points
      read(Input,*) NRPOINT
      write(*,*)'input_DEFAULT: READING POINTS...'
      write(*,*)'           NRPOINT = ',NRPOINT
!
      if (MAXNP.lt.NRPOINT) then
        write(*,7001) MAXNP,NRPOINT
 7001   format('input_POINTS: MAXNP,NRPOINT = ',2i6)
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'input_POINTS: READING POINTS...'
      endif
!
!  ...read in and store the points data
      do np=1,NRPOINT
        read(Input,*)type
        POINTS(np)%Type  = type
!
        select case(type)
!
!  .....read in and store the coordinates for regular points
        case('Regular')
          allocate(POINTS(np)%Rdata(3))
          read(Input,*) POINTS(np)%Rdata(1:3)
!
        case('CoorNrm')
          allocate(POINTS(np)%Rdata(6))
          read(Input,*) POINTS(np)%Rdata(1:6)
!
!  .....read in and store the intersecting surfaces data and
!       starting point for implicit  points
        case('Implicit')
          allocate(POINTS(np)%Idata(3))
          read(Input,*) POINTS(np)%Idata(1:3)
          allocate(POINTS(np)%Rdata(3))
          read(Input,*) POINTS(np)%Rdata(1:3)
!
        case default
           write(*,*)'input_POINTS: WRONG POINT TYPE'
        end select
!
      enddo
!
!
      end subroutine input_POINTS
!
!
!
c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : this is a LEGACY routine. Do not use it, unless
c                    you revise it!
c//////////////////////////////////////////////////////////////////////
c
c
c----------------------------------------------------------------------
c
c   routine name       - def_impl
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 05
c
c   purpose            - routine defines implicit curves and
c                        rectangles based on the surface conformity
c                        information for hexahedrals
c
c
c   arguments
c     in:
c         Hexa_to_surf - surfaces to which hexa faces have to conform
c
c---------------------------------------------------------------------
c
      subroutine def_impl(Hexa_to_surf)
c
      use GMP
      use control
#include "syscom.blk"
c
      integer Hexa_to_surf(6,MAXHE)
c
c
c  ...for i-th edge,
c     edge_faces(1:4,i) contains face numbers defining the edge
      integer edge_faces(4,12)
c
c  ...for i-th edge
c     edge_to_vert(1:2,i) are vertex end points
      integer edge_to_vert(2,12)
c
c  ...for i-th edge
c     edge_to_face_vert(1:4,1:4,i) lists points on the four faces
c     defining the planes for the implicit definition of the edge
      integer edge_to_face_vert(4,4,12)
c
c  ...for i-th face
c     face_to_vert(1:4,i) are the face vertex numbers
      integer face_to_vert(4,6)
      dimension xnod(3,8)
c
c  ...auxiliary vectors
      dimension aux1(3),aux2(3),aux3(3),aux4(3)
c
c  ...four points on a face
      dimension xface_points(3,4)
c
      dimension void(3)
c
      data edge_faces /1,3,6,4, 1,4,3,5, 1,5,6,4, 1,6,3,5,
     .                 2,3,6,4, 2,4,3,5, 2,5,6,4, 2,6,3,5,
     .                 6,3,1,2, 3,4,1,2, 4,5,1,2, 5,6,1,2/
      data edge_to_vert /1,2, 2,3, 4,3, 1,4, 5,6, 6,7, 8,7, 5,8,
     .                   1,5, 2,6, 3,7, 4,8/
      data face_to_vert /1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6,
     .                   4,3,7,8, 1,4,8,5/
      data edge_to_face_vert /1,2,3,4, 1,2,6,5, 1,4,8,5, 2,3,7,6,
     .                        2,3,4,1, 2,3,7,6, 2,6,5,1, 3,7,8,4,
     .                        3,4,1,2, 4,3,7,8, 4,8,5,1, 3,7,6,2,
     .                        4,1,2,3, 1,4,8,5, 1,2,6,5, 4,3,7,8,
     .                        5,6,7,8, 6,5,1,2, 5,1,4,8, 6,2,3,7,
     .                        6,7,8,5, 7,6,2,3, 6,5,1,2, 7,8,4,3,
     .                        7,8,5,6, 7,8,4,3, 8,5,1,4, 7,6,2,3,
     .                        8,5,6,7, 8,5,1,4, 5,1,2,6, 8,4,3,7,
     .                        5,1,4,8, 5,1,2,6, 1,2,3,4, 5,6,7,8,
     .                        2,6,5,1, 6,2,3,7, 2,3,4,1, 6,7,8,5,
     .                        3,7,6,2, 3,7,8,4, 3,4,1,2, 7,8,5,6,
     .                        8,4,3,7, 4,8,5,1, 4,1,2,3, 8,5,6,7/
c
c---------------------------------------------------------------------
c
      iprint=0
c
c  ...loop through hexas
      do nh=1,NRHEXAS
        iflag=0
        do if=1,6
          if (Hexa_to_surf(if,nh).ne.0) iflag=1
        enddo
        if (iflag.eq.0) cycle
        HEXAS(nh)%Type = 'TrIntHex'
c
c  .....compute the vertex coordinates
        do jv=1,8
          np = HEXAS(nh)%VertNo(jv)
          call pointr(np, xnod(1:3,jv))
          if (iprint.eq.1) then
            write(*,7045) jv, xnod(1:3,jv)
 7045       format('def_impl: jv, xnod = ',i2,2x,3f14.6)
          endif
        enddo
c
c  .....check consistency
        do if=1,6
          nf = Hexa_to_surf(if,nh)
          if (nf.ne.0) then
            do ivf=1,4
              iv = face_to_vert(ivf,if)
              if (SURFACES(nf)%Type.eq.'Sphere') then
c
c  ..............a quick correction to patch up Adam's file
                 aux1(1:3) = xnod(1:3,iv)-SURFACES(nf)%Rdata(1:3)
                 s = sqrt(aux1(1)**2+aux1(2)**2+aux1(3)**2)
                 np = HEXAS(nh)%VertNo(iv)
                 POINTS(np)%Rdata(1:3) = SURFACES(nf)%Rdata(1:3) +
     .                aux1(1:3)/s*SURFACES(nf)%Rdata(4)
                 xnod(1:3,iv) = POINTS(np)%Rdata(1:3)
              endif
c
              call surf(nf,xnod(1:3,iv), fval,void)
              if (abs(fval).gt.GEOM_TOL) then
                write(*,8001) nh,if,ivf,iv,xnod(1:3,iv),nf,
     .                        SURFACES(nf)%Type,fval,GEOM_TOL
 8001           format('def_impl: INCONSISTENT DATA ',
     .               /,' hexa   NUMBER = ',i5,
     .               /,' face   NUMBER = ',i5,
     .               /,' face vertex and hexa vertex NUMBER = ',2i5,
     .               /,' point COORDINATES = ',3e12.5,
     .               /,' surface NUMBER = ',i5,
     .               /,' surface TYPE = ',a,
     .               /,' fval = ',e12.5,
     .               /,' GEOM_TOL = ',e12.5)
                stop 1
              endif
            enddo
          endif
        enddo
c
c  .....loop through edges
        do ie=1,12
          if1 = edge_faces(1,ie); if2 = edge_faces(2,ie)
          nf1 = Hexa_to_surf(if1,nh); nf2 = Hexa_to_surf(if2,nh)
          if3 = edge_faces(3,ie); if4 = edge_faces(4,ie)
          nf3 = Hexa_to_surf(if3,nh); nf4 = Hexa_to_surf(if4,nh)
c
          if ((nf1.ne.0).or.(nf2.ne.0)) then
c
c  .........if the curve is still defined as a segment of line
c           and the orientation is consistent with the local
c           one (it must be for at least one adjacent block
c           that was used to defined the curve)...
            nc = HEXAS(nh)%EdgeNo(ie)
c
c  .........endpoints vertex numbers
            if (nc.gt.0) then
              iv1 = edge_to_vert(1,ie); iv2 = edge_to_vert(2,ie)
            else
              iv1 = edge_to_vert(2,ie); iv2 = edge_to_vert(1,ie)
            endif
            nc = abs(nc)
            if (CURVES(nc)%Type.ne.'Seglin') cycle
            if (iprint.eq.1) then
              write(*,7004) nh,ie
 7004         format('def_impl: REDEFINING CURVE,'
     .               ' HEXA, EDGE = ',i6,i3)
            endif
c
c  .........a curve lying on a single surface
            if ((nf1.eq.0).or.(nf2.eq.0)) then
              nf = nf1+nf2
              select case (SURFACES(nf)%Type)
c
              case('Cylinder')
c
c  .............proceed to define a cylinder geodesics
                CURVES(nc)%Type  = 'CylGeod'
                allocate(CURVES(nc)%Idata(1), STAT = ic)
                CURVES(nc)%Idata(1) = nf
                allocate(CURVES(nc)%Rdata(12), STAT = ic)
c
c  .............generate the data for an explicit parametrization of the
c               cylinder geodesics
                call gen_Cyl_Geod(nf,xnod(1:3,iv1),xnod(1:3,iv2),
     .          CURVES(nc)%Rdata(1),
     .          CURVES(nc)%Rdata(2),CURVES(nc)%Rdata(3),
     .          CURVES(nc)%Rdata(4:12))
                if (iprint.eq.1) then
                  write(*,7001) nc, CURVES(nc)%Idata(1)
 7001            format('def_impl: HAVE DEFINED CYLINDER GEODESICS ',i4,
     .                   ' CYLINDER = ',i4)
                  call pause
                endif
c
              case('Sphere')
c
c  .............save the sphere number
                nff = nf
c
c  .............use implict curve to define geodesics
                CURVES(nc)%Type  = 'ImpCir'
                allocate(CURVES(nc)%Idata(4), STAT = ic)
c
c  .............first surface is the sphere
                CURVES(nc)%Idata(1) = nff

c
c  .............second surface is a plane passing through the two endpoints
c               and the center of the sphere
                NRSURFS = NRSURFS + 1
                if (NRSURFS.gt.MAXSU) then
                  write(*,7003) NRSURFS,MAXSU
 7003             format('def_impl: MAXSU EXCEEDED, NRSURFS,MAXSU = ',
     .                    2i7)
                  stop 1
                endif
                nf = NRSURFS
                SURFACES(nf)%Type = 'ThrPt'
                allocate(SURFACES(nf)%Rdata(9), STAT = ic)
                SURFACES(nf)%Rdata(1:3) = SURFACES(nff)%Rdata(1:3)
                SURFACES(nf)%Rdata(4:6) = xnod(1:3,iv1)
                SURFACES(nf)%Rdata(7:9) = xnod(1:3,iv2)
                CURVES(nc)%Idata(2) = nf
c
c  .............third surface is a plane perpendicular to the tangent
c               of the great cicrle at first endpoint
                aux1(1:3) = xnod(1:3,iv1) - SURFACES(nff)%Rdata(1:3)
                aux2(1:3) = xnod(1:3,iv2) - xnod(1:3,iv1)
                aux3(1) = aux1(2)*aux2(3) - aux1(3)*aux2(2)
                aux3(2) = aux1(3)*aux2(1) - aux1(1)*aux2(3)
                aux3(3) = aux1(1)*aux2(2) - aux1(2)*aux2(1)
c
                aux4(1) = aux3(2)*aux1(3) - aux3(3)*aux1(2)
                aux4(2) = aux3(3)*aux1(1) - aux3(1)*aux1(3)
                aux4(3) = aux3(1)*aux1(2) - aux3(2)*aux1(1)
                s = sqrt(aux4(1)**2+aux4(2)**2+aux4(3)**2)
                NRSURFS = NRSURFS + 1
                if (NRSURFS.gt.MAXSU) then
                  write(*,7003) NRSURFS,MAXSU
                  stop 1
                endif
                nf = NRSURFS
                SURFACES(nf)%Type = 'VecPt'
                allocate(SURFACES(nf)%Rdata(6), STAT = ic)
                SURFACES(nf)%Rdata(1:3) = xnod(1:3,iv1)
                SURFACES(nf)%Rdata(4:6) = aux4(1:3)/s
                CURVES(nc)%Idata(3) = nf
c
c  .............same for the fourth surface but at the other endpoint
                aux1(1:3) = xnod(1:3,iv2) - SURFACES(nff)%Rdata(1:3)
                aux2(1:3) = xnod(1:3,iv2) - xnod(1:3,iv1)
                aux3(1) = aux1(2)*aux2(3) - aux1(3)*aux2(2)
                aux3(2) = aux1(3)*aux2(1) - aux1(1)*aux2(3)
                aux3(3) = aux1(1)*aux2(2) - aux1(2)*aux2(1)
c
                aux4(1) = aux3(2)*aux1(3) - aux3(3)*aux1(2)
                aux4(2) = aux3(3)*aux1(1) - aux3(1)*aux1(3)
                aux4(3) = aux3(1)*aux1(2) - aux3(2)*aux1(1)
                s = sqrt(aux4(1)**2+aux4(2)**2+aux4(3)**2)
                NRSURFS = NRSURFS + 1
                if (NRSURFS.gt.MAXSU) then
                  write(*,7003) NRSURFS,MAXSU
                  stop 1
                endif
                nf = NRSURFS
                SURFACES(nf)%Type = 'VecPt'
                allocate(SURFACES(nf)%Rdata(6), STAT = ic)
                SURFACES(nf)%Rdata(1:3) = xnod(1:3,iv2)
                SURFACES(nf)%Rdata(4:6) = aux4(1:3)/s
                CURVES(nc)%Idata(4) = nf
                if (iprint.eq.1) then
                  write(*,7005) nc, CURVES(nc)%Idata(1:4)
 7005             format('def_impl: HAVE DEFINED SS IMPLICIT CURVE', i4,
     .                   ' SURFACES = ',4i4)
c
c  ...............check compatibility
                  do i=1,4
                    nf = CURVES(nc)%Idata(i)
                    if (i.ne.4) then
                      call surf(nf,xnod(1:3,iv1), fval,aux1)
                      if (abs(fval).gt.GEOM_TOL) then
                        write(*,*) 'iv,i,nf,fval = ',iv1,i,nf,fval
                        stop 1
                      endif
                    endif
                    if (i.ne.3) then
                      call surf(nf,xnod(1:3,iv2), fval,aux1)
                      if (abs(fval).gt.GEOM_TOL) then
                        write(*,*) 'iv,i,nf,fval = ',iv2,i,nf,fval
                        stop 1
                      endif
                    endif
                  enddo
                  call pause
                endif
c
              case('Cone')
c
c  .............save the cone number
                nff = nf
c
c  .............use implict curve to define the curve
                CURVES(nc)%Type  = 'ImpCir'
                allocate(CURVES(nc)%Idata(4), STAT = ic)
c
c  .............first surface is the cone
                CURVES(nc)%Idata(1) = nff
c
c  .............second surface is a plane passing through the two endpoints
c               and the projection of the mid-point onto the cone axis

                NRSURFS = NRSURFS + 1
                if (NRSURFS.gt.MAXSU) then
                  write(*,7003) NRSURFS,MAXSU
                  stop 1
                endif
                nf = NRSURFS
                SURFACES(nf)%Type = 'ThrPt'
                allocate(SURFACES(nf)%Rdata(9), STAT = ic)
c
c  .............mid-point of the segment
                aux1(1:3) = (xnod(1:3,iv1)+xnod(1:3,iv2))/2.d0
c
c  .............position vector of the mid-point
                aux1(1:3) = aux1(1:3) - SURFACES(nff)%Rdata(1:3)
c
c  .............unit vector of the cone vector
                s = 0.d0
                do jj=1,3
                  s = s + SURFACES(nff)%Rdata(3+jj)**2
                enddo
                s = sqrt(s)
                aux2(1:3) = SURFACES(nff)%Rdata(4:6)/s
c
c  .............cone zeta coordinate of the mid-point
                zet = aux1(1)*aux2(1)+aux1(2)*aux2(2)+aux1(3)*aux2(3)
c
c  .............projection
                aux3(1:3) = SURFACES(nff)%Rdata(1:3)
     .                    + zet*aux2(1:3)
c
c  .............define the plane
                SURFACES(nf)%Rdata(1:3) = aux3(1:3)
                SURFACES(nf)%Rdata(4:6) = xnod(1:3,iv1)
                SURFACES(nf)%Rdata(7:9) = xnod(1:3,iv2)
                CURVES(nc)%Idata(2) = nf
c
c  .............third surface is a plane perpendicular to the segment
                aux1(1:3) = xnod(1:3,iv2)-xnod(1:3,iv1)
                s = sqrt(aux1(1)**2+aux1(2)**2+aux1(3)**2)
                NRSURFS = NRSURFS + 1
                if (NRSURFS.gt.MAXSU) then
                  write(*,7003) NRSURFS,MAXSU
                  stop 1
                endif
                nf = NRSURFS
                SURFACES(nf)%Type = 'VecPt'
                allocate(SURFACES(nf)%Rdata(6), STAT = ic)
                SURFACES(nf)%Rdata(1:3) = xnod(1:3,iv1)
                SURFACES(nf)%Rdata(4:6) = aux1(1:3)/s
                CURVES(nc)%Idata(3) = nf
c
c  .............same for the fourth surface but at the other endpoint
                NRSURFS = NRSURFS + 1
                if (NRSURFS.gt.MAXSU) then
                  write(*,7003) NRSURFS,MAXSU
                  stop 1
                endif
                nf = NRSURFS
                SURFACES(nf)%Type = 'VecPt'
                allocate(SURFACES(nf)%Rdata(6), STAT = ic)
                SURFACES(nf)%Rdata(1:3) = xnod(1:3,iv2)
                SURFACES(nf)%Rdata(4:6) = aux1(1:3)/s
                CURVES(nc)%Idata(4) = nf
                if (iprint.eq.1) then
                  write(*,7015) nc, CURVES(nc)%Idata(1:4)
 7015             format('def_impl: HAVE DEFINED SC IMPLICIT CURVE', i4,
     .                   ' SURFACES = ',4i4)
c
c  ...............check compatibility
                  do i=1,4
                    nf = CURVES(nc)%Idata(i)
                    if (i.ne.4) then
                      call surf(nf,xnod(1:3,iv1), fval,aux1)
                      if (abs(fval).gt.GEOM_TOL) then
                        write(*,*) 'iv,i,nf,fval = ',iv1,i,nf,fval
                        stop 1
                      endif
                    endif
                    if (i.ne.3) then
                      call surf(nf,xnod(1:3,iv2), fval,aux1)
                      if (abs(fval).gt.GEOM_TOL) then
                        write(*,*) 'iv,i,nf,fval = ',iv2,i,nf,fval
                        stop 1
                      endif
                    endif
                  enddo
                  call pause
                endif
              case default
                write(*,*) 'def_impl: UNFINISHED 1'
                stop 1
              end select
c
c  .........intersection of two curvilinear surfaces
            else
c
c  ...........use the implicit curve, first two surfaces are known
              CURVES(nc)%Type  = 'ImpCir'
              allocate(CURVES(nc)%Idata(4), STAT = ic)
              CURVES(nc)%Idata(1) = nf1
              CURVES(nc)%Idata(2) = nf2
c
c  ...........define the two remaining surfaces as planes
              NRSURFS = NRSURFS + 1
              if (NRSURFS.gt.MAXSU) then
                write(*,7003) NRSURFS,MAXSU
                stop 1
              endif
              nf = NRSURFS
              SURFACES(nf)%Type = 'VecPt'
              allocate(SURFACES(nf)%Rdata(6), STAT = ic)
              SURFACES(nf)%Rdata(1:3) = xnod(1:3,iv1)
              do jj=1,4
                iv = edge_to_face_vert(jj,3,ie)
                xface_points(1:3,jj) = xnod(1:3,iv)
              enddo
              call determine_plane1(xnod(1:3,iv1),
     .             xface_points(1:3,2:4),SURFACES(nf)%Rdata(4:6))
              CURVES(nc)%Idata(3) = nf
              if (iprint.eq.1) then
                write(*,7006) nf,SURFACES(nf)%Rdata(1:6)
 7006           format('def_impl: HAVE DEFINED PLANE, nf = ',i5,
     .                 ' POINT = ',3e12.5,' RN = ',3e12.5)
              endif
c
              NRSURFS = NRSURFS + 1
              if (NRSURFS.gt.MAXSU) then
                write(*,7003) NRSURFS,MAXSU
                stop 1
              endif
              nf = NRSURFS
              SURFACES(nf)%Type = 'VecPt'
              allocate(SURFACES(nf)%Rdata(6), STAT = ic)
              SURFACES(nf)%Rdata(1:3) = xnod(1:3,iv2)
              do jj=1,4
                iv = edge_to_face_vert(jj,4,ie)
                xface_points(1:3,jj) = xnod(1:3,iv)
              enddo
              call determine_plane1(xnod(1:3,iv2),
     .             xface_points(1:3,2:4),SURFACES(nf)%Rdata(4:6))
              CURVES(nc)%Idata(4) = nf
              if (iprint.eq.1) then
                write(*,7006) nf,SURFACES(nf)%Rdata(1:6)
              endif
              if (iprint.eq.1) then
                write(*,7012) nc, CURVES(nc)%Idata(1:4)
 7012           format('def_impl: HAVE DEFINED DS IMPLICIT CURVE', i4,
     .                   ' SURFACES = ',4i4)
                call pause
c
c  .............check compatibility
                do i=1,4
                  nf = CURVES(nc)%Idata(i)
                  if (i.ne.4) then
                    call surf(nf,xnod(1:3,iv1), fval,aux1)
                    if (abs(fval).gt.GEOM_TOL) then
                      write(*,*) 'iv,i,nf,fval = ',iv1,i,nf,fval
                      stop 1
                    endif
                  endif
                  if (i.ne.3) then
                    call surf(nf,xnod(1:3,iv2), fval,aux1)
                    if (abs(fval).gt.GEOM_TOL) then
                      write(*,*) 'iv,i,nf,fval = ',iv2,i,nf,fval
                      stop 1
                    endif
                  endif
                enddo
                call pause
              endif
c
            endif
c
c  .........modify the definition of all adjacent rectas
            do if=1,CURVES(nc)%NrFig
              call decode(abs(CURVES(nc)%FigNo(if)), nr,lab)
              select case(lab)
              case(2)
                RECTANGLES(nr)%Type = 'TraQua'
              case default
                write(*,*) 'def_impl: UNFINISHED 2'
                stop 1
              end select
            enddo
c
c  .......end of implicit curve definition
          endif
c
c  .....end of loop through edges
        enddo
c
c  ...end of loop through hexahedra
      enddo
c
c---------------------------------------------------------------------
c
c  ...loop through hexas
      do nh=1,NRHEXAS
c
c  .....loop through faces
        do jf=1,6
c
c  .......if the face has to conform to a surface
          ns = hexa_to_surf(jf,nh)
          if (ns.ne.0) then
            call decode(HEXAS(nh)%FigNo(jf), nr,no)
c
c  .........if the orientation is 0 (it must be at least for one
c           adjacent hexa! )
            if (no.eq.0) then
              if (iprint.eq.2) then
                write(*,7025) nh,jf,ns,nr,
     .          (CURVES(abs(RECTANGLES(nr)%EdgeNo(ie)))%type, ie=1,4)
 7025           format('def_impl: nh,jf,ns,nr, CURVES TYPES = ',
     .                  i4,2x,i2,i5,i5,2x,4a10)
                call pause
              endif
c
c  ...........decide on the type of rectangle
              idec=1
c
c  ...........loop through the edge curves of the rectangle
              do ie=1,4
                nc = abs(RECTANGLES(nr)%EdgeNo(ie))
c
c  .............check if the curve is implicit
                if (CURVES(nc)%Type.ne.'ImpCir') then
                  idec=0
                  exit
                endif
c
c  ...........end of loop through the edge curves
              enddo
c
              select case(idec)
c
c  ...........conforming rectangle
              case(0)
                RECTANGLES(nr)%Type = 'ConfRec'
                allocate(RECTANGLES(nr)%Idata(1), STAT = ic)
                RECTANGLES(nr)%Idata(1) = ns
c
c  ...........implicit rectangle
              case(1)
                RECTANGLES(nr)%Type = 'ImpRec'
                allocate(RECTANGLES(nr)%Idata(5), STAT = ic)
                RECTANGLES(nr)%Idata(1) = ns
                do ie=1,4
                  nc = abs(RECTANGLES(nr)%EdgeNo(ie))
c
c  ...............recycle surface used to defined the implicit curve
c                 to define the implicit rectangle..................
c
c  ...............get the mid-point of the edge curve
                  xi = 0.5d0
                  call curve(nc,xi, aux4,aux1)
c
c  ...............compute gradients for the main surface defining
c                 the rectangle and the two surfaces defining the curve
                  call surf(ns,aux4, fval3,aux3)
                  call surf(CURVES(nc)%Idata(1),aux4, fval1,aux1)
                  call surf(CURVES(nc)%Idata(2),aux4, fval2,aux2)
                  call cross_product(aux3,aux1, aux4)
                  fval1 = aux4(1)**2+aux4(2)**2+aux4(3)**2
                  call cross_product(aux3,aux2, aux4)
                  fval2 = aux4(1)**2+aux4(2)**2+aux4(3)**2
                  if (fval1.gt.fval2) then
                    RECTANGLES(nr)%Idata(1+ie) = CURVES(nc)%Idata(1)
                  else
                    RECTANGLES(nr)%Idata(1+ie) = CURVES(nc)%Idata(2)
                  endif
c
cldem Dec 15 05
ccc                  if (CURVES(nc)%Idata(1).eq.ns) then
ccc                    RECTANGLES(nr)%Idata(1+ie) = CURVES(nc)%Idata(2)
ccc                  else
ccc                    RECTANGLES(nr)%Idata(1+ie) = CURVES(nc)%Idata(1)
ccc                  endif
c
c  .............end of loop through the edge curves
                enddo
                if (iprint.eq.1) then
                  write(*,7002) nr, RECTANGLES(nr)%Idata(1:5)
 7002             format('def_impl: DEFINED IMPLICIT RECTA ', i4,
     .                   ' SURFACES = ',5i4)
                  call pause
                endif
c
              end select
c
c  .........if the right orientation
            endif
c
c  .......if the face has to conform to a surface
          endif
c
c  .....end of loop through faces
        enddo
c
c  ...end of loop through hexas
      enddo
c
c
      end subroutine def_impl
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - gen_Cyl_Geod
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 05
c
c   purpose            - given a cylinder and two points on it,
c                        routine evaluates data necessary for
c                        determining a geodesics connecting the two points
c
c   arguments :
c     in:
c               Ns     - a surface number (must be a cylinder)
c               X1,X2  - coordinates of two points lying on the
c                        cylinder
c     out:
c               Theta2 - local theta coordinate for the endpoint
c                        of the curve (Theta1 = 0)
c               Rz1,Rz2  - z coordinates of the endpoints of the curve
c               Transf - transformation matrix from a local,
c                        cylinder-system to the global Cartesian
c                        system of coordinates
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine gen_Cyl_Geod(Ns,X1,X2, Theta2,Rz1,Rz2,Transf)
c
      use GMP
      use control
#include "syscom.blk"
c
      dimension X1(3),X2(3), Transf(3,3)
c
c  ...cylinder origin and axis, auxiliary space
      dimension xcenter(3),xaxis(3),aux(3)
c
c  ...endpoints coordinates in the local Cartesian system
      dimension x1p(3),x2p(3)
c
      data epsilon /1.d-13/
c
      iprint=0
c
      pihalf = dacos(0.d0)
c
      if (iprint.eq.1) then
        write(*,7000) X1,X2
 7000   format('gen_Cyl_Geod: X1,X2 = ',2(3e12.5,2x))
      endif
c
c  ...check the type
      if (SURFACES(Ns)%Type.ne.'Cylinder') then
        write(*,7001) ns,SURFACES(ns)%Type
 7001   format('gen_Cyl_Geod: WRONG SURFACE TYPE, No, Type = ', i4,a10)
        stop 1
      endif
      xcenter(1:3) = SURFACES(ns)%Rdata(1:3)
      xaxis(1:3)   = SURFACES(ns)%Rdata(4:6)
      radius   = SURFACES(ns)%Rdata(7)
c
c  ...check if the two points are indeed on the cylinder
      call surf(Ns,X1, fval,aux)
      if (abs(fval).gt.GEOM_TOL) then
        write(*,7002) X1,Ns,fval
 7002   format('gem_Cyl_Geod: POINT ',3e12.5,' NOT ON SURFACE ',i4,
     .         'fval = ',e12.5)
        stop 1
      endif
      call surf(Ns,X2, fval,aux)
      if (abs(fval).gt.GEOM_TOL) then
        write(*,7002) X2,Ns,fval
        stop 1
      endif
c
c  ...the third local cylinder versor coincides with its axis
      s = sqrt(xaxis(1)**2+xaxis(2)**2+xaxis(3)**2)
      Transf(1:3,3) = xaxis(1:3)/s
c
c  ...compute cross product of the cylinder vector with the vector
c     connecting its origin with the first curve endpoint
      aux(1:3) = X1(1:3) - xcenter(1:3)
      call cross_product(xaxis,aux, Transf(1:3,2))
c
c  ...normalize
      s = sqrt(Transf(1,2)**2+Transf(2,2)**2+Transf(3,2)**2)
      Transf(1:3,2) = Transf(1:3,2)/s
c
c  ...compute the first axis as the cross product of the second and third
c     axes versors
      call cross_product(Transf(1:3,2),Transf(1:3,3), Transf(1:3,1))
c
      if (iprint.eq.1) then
        do i=1,3
          write(*,7003) i,Transf(1:3,i)
 7003     format('gen_Cyl_Geod: ',i1,'-TH VERSOR = ',3f8.3)
        enddo
      endif
c
c  ...compute local Cartesian coordinates
      call transform(12,Transf,aux,x1p)
      aux(1:3) = X2(1:3) - xcenter(1:3)
      call transform(12,Transf,aux,x2p)
c
      if (iprint.eq.1) then
        write(*,7004) x1p,x2p
 7004   format('gen_Cyl_Geod: x1p,x2p = ',2(3e12.5,2x))
      endif
c
      Rz1 = x1p(3); Rz2 = x2p(3)
      if (abs(x1p(2)).gt.GEOM_TOL) then
        write(*,7005) x1p
 7005   format('gen_Cyl_Geod: INCONSISTENCY ! x1p = ',3e12.5)
        stop 1
      endif
      call cart_to_polar(x2p, void,Theta2)
c
      if (iprint.eq.1) then
        write(*,7006) Theta2,Rz1,Rz2
 7006   format('gen_Cyl_Geod: Theta2,Rz1,Rz2 = ',e12.5,2x,2e12.5)
        call pause
      endif
c
c
      end subroutine gen_Cyl_Geod
