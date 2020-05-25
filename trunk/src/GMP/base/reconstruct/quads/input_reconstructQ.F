c----------------------------------------------------------------------
c
c   routine name       - input_reconstruct
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Nov 07
c
c   purpose            - Routine reads in the following data,
c                          number of points on the surface
c                          and number of points total
c                            for each point on the surface:
c                              coordinates of the point and
c                              components of a normal to the surface
c                              at the point
c                            for each point in the volume (including
c                            those on the surface):
c                              coordinates of the point
c                          number of hexas
c                            for each hexa element in the volume:
c                              element-to-point connectivities
c                        The routine generates then a full GMP
c                        information, defining curves and remaining
c                        rectangles and constructing the full GMP
c                        connectivities. The rectangles on the surface
c                        undergo the geometry reconstruction discussed
c                        in the book, and the hexas use transfinite
c                        interpolation
c
c   arguments            none
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine input_reconstruct
c
      use control
      use GMP
#include "syscom.blk"
#include "cinout.blk"
c
      integer, dimension(:), pointer:: vpoint_to_spoint
      integer, dimension(:,:), pointer:: point_to_blocks
      integer, dimension(:),   pointer:: point_nrbl
c
c  ...work space
      dimension rdata(6)
c
      parameter (mphex = 20)
c
      integer   hexa_edge_to_vert(2,12),
     .          hexa_face_to_vert(4,6),
     .          hexa_vert_to_edge(3,8),
     .          hexa_vert_to_face(3,8),
     .          list(2,20)
c     for i-th block on the list
c       list(1,i) = block connected to an edge
c       list(2,i) = the corresponding edge number
c     or
c       list(1,i) = block connected to a face
c       list(2,i) = the corresponding face number
c
      character(10) type
c
      data hexa_edge_to_vert /1,2, 2,3, 4,3, 1,4, 5,6, 6,7, 8,7, 5,8,
     .                        1,5, 2,6, 3,7, 4,8/
      data hexa_face_to_vert /1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6,
     .                        4,3,7,8, 1,4,8,5/
      data hexa_vert_to_edge /4,1,9, 1,2,10, 2,3,11, 3,4,12,
     .                        8,5,9, 5,6,10, 6,7,11, 7,8,12/
      data hexa_vert_to_face /1,6,3, 1,3,4, 1,4,5, 1,5,6,
     .                        2,6,3, 2,3,4, 2,4,5, 2,5,6/
c
c---------------------------------------------------------------------
c
      iprint=10
c
c  ...set up the dimension of the problem
      NDIM=3
c
c  ...set up the unused parameters
      NRSURFS=0; NRTRIAN=0; NRPRISM=0
c
c---------------------------------------------------------------------
c
c  Step 1: Input points
c
c  ...read the number of points on the surface and in the volume
      read(KIN,*) nrsp, NRPOINT
c
c  ...allocate the memory dynamically
      allocate(POINTS(NRPOINT), STAT=is)
      if (is.ne.0) then
        write(*,*) 'input_reconstruct: POINTS ARE NOT ALLOCATED'
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'input_reconstruct: READING SURFACE POINTS...'
      endif
c
c  ...read in the surface points
      do np=1,nrsp
        POINTS(np)%Type  = 'CoorNrm'
        POINTS(np)%NrCurv = 0
        allocate(POINTS(np)%Rdata(6), STAT=is)
        if (is.ne.0) then
          write(*,*) 'input_reconstruct: Rdata IS  NOT ALLOCATED',
     .               'POINT NUMBER = ', np
          stop 1
        endif
        read(KIN,*) POINTS(np)%Rdata(1:6)
c
c  .....renormalize the normal
        call normalize(POINTS(np)%Rdata(4:6))
      enddo
c
c  ...allocate volume points to surface points connectivities
      allocate(vpoint_to_spoint(NRPOINT), STAT=is)
      if (is.ne.0) then
        write(*,*) 'input_reconstruct: vpoint_to_spoint ALLOCATION ',
     .             ' FAILED'
        stop 1
      endif
      vpoint_to_spoint(1:NRPOINT)=0
c
      if (iprint.eq.1) then
        write(*,*) 'input_reconstruct: READING VOLUME POINTS'
      endif
c
c  ...read in the volume points
      np=nrsp
      do ip=1,NRPOINT
        read(KIN,*) rdata(1:3)
c
c  .....search for a point on the surface
        do jp=1,nrsp
          dist=0.d0
          do ivar=1,3
            dist = dist + (POINTS(jp)%Rdata(ivar)-rdata(ivar))**2
          enddo
          dist = sqrt(dist)
c
c  .......have found a matching point on the surface
          if (dist.lt.GEOM_TOL) then
c
c  .........save the connectivity
            vpoint_to_spoint(ip)=jp
          endif
c
c  .....end of loop through surface points
        enddo
c
c  .....if no surface point has been found...
        if (vpoint_to_spoint(ip).eq.0) then
c
c  .......generate a new regular point
          np=np+1
          if (np.gt.NRPOINT) then
            write(*,8001) np,NRPOINT
 8001       format('input_reconstruct: INCONSISTENCY 1, np,NRPOINT = ',
     .             2i6)
            stop 1
          endif
          POINTS(np)%Type  = 'Regular'
          POINTS(np)%NrCurv = 0
          allocate(POINTS(np)%Rdata(3), STAT=is)
          if (is.ne.0) then
            write(*,*) 'input_reconstruct: Rdata IS  NOT ALLOCATED',
     .                 'POINT NUMBER = ', np
            stop 1
          endif
          POINTS(np)%Rdata(1:3) = rdata(1:3)
c
c  .......save the connectivity
          vpoint_to_spoint(ip)=np
        endif
c
c  ...end of loop through volume points
      enddo
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: HAVE READ POINTS'
      endif
c
c----------------------------------------------------------------------
c
c  Step 2: Input hexas
c
c  ...read the number of hexas in the volume
      read(KIN,*) NRHEXAS
      if (iprint.eq.1) then
        write(*,*) 'input_reconstruct: NRHEXAS = ',NRHEXAS
      endif
c
c  ...set up the dimension of the manifold
      if (NRHEXAS.le.0) then
        write(*,8002) NRHEXAS
 8002   format('input_reconstruct: NRHEXAS = ',i6)
        stop 1
      else
        MANDIM=3
      endif
c
c  ...allocate dynamically the memory
      allocate(HEXAS(NRHEXAS), STAT=is )
      if (is.ne.0) then
        write(*,*) 'input_reconstruct: HEXAS ARE NOT ALLOCATED'
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'input_reconstruct: READING HEXAS...'
      endif
c
      do ih=1,NRHEXAS
        HEXAS(ih)%Type = 'TriLiHex'
        read(KIN,*) HEXAS(ih)%VertNo(1:8)
c
c  .....change the enumeration
        do i=1,8
          j = HEXAS(ih)%VertNo(i)
          HEXAS(ih)%VertNo(i) = vpoint_to_spoint(j)
        enddo
        if (iprint.eq.1) then
          write(*,7636) ih, HEXAS(ih)%VertNo(1:8)
 7636     format('HEXA = ',i6,' VERTICES = ',8i7)
        endif
c
c  ...end of loop through hexahedra
      enddo
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: HAVE READ HEXAS'
        call pause
      endif
c
c----------------------------------------------------------------------
c
c  Step 3: Determine point to hexas connectivities
c
      allocate(point_to_blocks(1:mphex,1:NRPOINT), STAT=is )
      allocate(point_nrbl(1:NRPOINT), STAT=is)
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
            nrbl = nrbl+1
            if (nrbl.gt.mphex) then
              write(*,*) 'input_reconstruct: INCREASE mphex = ',mphex
              stop 1
            endif
            point_to_blocks(nrbl,np) = nick
            point_nrbl(np) = nrbl
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
 7005     format('input_reconstruct: np = ',i7,
     .           ' CONNECTED BLOCKS = ',15i8)
        enddo
        call pause
      endif
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: HAVE CONNECTED HEXAS TO POINTS'
      endif
c
c----------------------------------------------------------------------
c
c  Step 4a: Count curves and determine temporary connectivities
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
            ibl=1
            list(1,ibl) = nick; list(2,ibl) = je
c
c  .........determine the endpoints of the edge
            iv1 = hexa_edge_to_vert(1,je)
            np1 = HEXAS(nh)%VertNo(iv1)
            iv2 = hexa_edge_to_vert(2,je)
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
                  ie = hexa_vert_to_edge(l,kv1)
c
c  ...............look for the other vertex point
                  do ive=1,2
                    iv3 = hexa_edge_to_vert(ive,ie)
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
                write(*,*) 'input_reconstruct: WRONG lab = ',lab
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
 7003         format('input_reconstruct: HAVE DEFINED CURVE ',i4,
     .               ' FOR HEXA ', i4,' AND EDGE ',i2)
              write(*,*) 'CONNECTED BLOCKS WITH EDGE NUMBERS'
              write(*,7007) (list(1,jj)/10,jj=1,ibl)
 7007         format(10i5)
              write(*,7007) (list(2,jj),jj=1,ibl)
ccc              call pause
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
      allocate(CURVES(NRCURVE), STAT=is)
      if (is.ne.0) then
        write(*,*) 'input_reconstruct: CURVES ARE NOT ALLOCATED'
        stop 1
      endif
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
            if (iprint.eq.1) then
              write(*,*) 'input_reconstruct: DEFINING CURVE nc = ',nc
            endif
c
c  .........determine the endpoints of the edge
            iv1 = hexa_edge_to_vert(1,je)
            np1 = HEXAS(nh)%VertNo(iv1)
            iv2 = hexa_edge_to_vert(2,je)
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
 7001     format('input_reconstruct: nc,points = ',i4,2x,2i5)
        enddo
        call pause
      endif
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: HAVE DEFINED AND CONNECTED',
     .             ' CURVES'
      endif
c
c----------------------------------------------------------------------
c
c  Step 5: count rectas and determine temporary hexas to face rectas
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
            iv = hexa_face_to_vert(1,jf)
            np1 = HEXAS(nh)%VertNo(iv)
            iv = hexa_face_to_vert(2,jf)
            np2 = HEXAS(nh)%VertNo(iv)
            iv = hexa_face_to_vert(3,jf)
            np3 = HEXAS(nh)%VertNo(iv)
            iv = hexa_face_to_vert(4,jf)
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
                  if = hexa_vert_to_face(l,iv1)
c
c  ...............look for two two other common points
                  iflag=0
                  do ivf=1,4
                    iv = hexa_face_to_vert(ivf,if)
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
                write(*,*) 'input_reconstruct: WRONG lab = ',lab
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
            if (iprint.eq.5) then
              write(*,7008) nr,nh,jf
 7008         format('input_reconstruct: HAVE DEFINED RECTA ',i4,
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
      allocate( RECTANGLES(NRRECTA), STAT=is )
      if (is.ne.0) then
        write(*,*) 'input_reconstruct: RECTANGLES ARE NOT ALLOCATED'
        stop 1
      endif
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
            iv = hexa_face_to_vert(1,jf)
            np1 = HEXAS(nh)%VertNo(iv)
            iv = hexa_face_to_vert(2,jf)
            np2 = HEXAS(nh)%VertNo(iv)
            iv = hexa_face_to_vert(3,jf)
            np3 = HEXAS(nh)%VertNo(iv)
            iv = hexa_face_to_vert(4,jf)
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
 7002     format('input_reconstruct: nr,points = ',i4,2x,4i5)
        enddo
        call pause
      endif
c
c----------------------------------------------------------------------
c
c  Step 6: complete connectivities
c
      deallocate(point_to_blocks,point_nrbl, STAT=is)
      deallocate(vpoint_to_spoint, STAT=is)
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: calling connect'
      endif
      call connect
c
c  ...data consistency check
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: CHECKING..'
      endif
ccc      call graphg
ccc      call check_reconstruct
c
c
c----------------------------------------------------------------------
c
c  Step 7: reconstruct the geometry
c
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: BEGINNING RECONSTRUCTION'
      endif
      call reconstruct
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: FINISHED RECONSTRUCTION'
      endif
c
c
      end

c----------------------------------------------------------------------
c
c   routine name       - check_reconstruct
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - May 07
c
c   purpose            - Routine verifies orientation and determinants
c                        for a hexahedral mesh with some points on
c                        an oriented surface
c
c   arguments            none
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine check_reconstruct
c
      use control
      use GMP
#include "syscom.blk"
#include "cinout.blk"
c
c  ...vertex coordinates, normals, surface flags
      dimension xv(3,8),rn(3,8),nflag(8)
c
c  ...normal for a face vertex
      dimension rnf(3)
c
c  ...order of approximation (linear)
      dimension norder(19)
c
c  ...shape functions and their derivatives
      dimension vshap(8),dvshap(3,8)
c
c  ...master element coordinates of vertices
      dimension xi(1:3,8)
c
c  ...derivatives wrt reference coordinates
      dimension dxdxi(3,3,8)
c
c  ...vertex points on a face, face coordinates
      dimension nface_vert(4,6), nface_axes(2,6)
c
c  ...neigboring vertices for a vertex in the cube
      dimension neig_vert(3,8)
c
c  ...a flag indicating whether a point has been visited
      integer, dimension(:), pointer:: nvisit
c
c  ...a work space
      dimension temp(3),temp1(3),iface(6),
     .          aa(3,3),bb(3),uu(3),ipiv(3)
c
      data norder /1,1,1,1,1,1,1,1,1,1,1,1,
     .             11,11,11,11,11,11, 111/
c
      data neig_vert/2,4,5, 3,1,6, 4,2,7, 1,3,8,
     .               8,6,1, 5,7,2, 6,8,3, 7,5,4/
c
      data nface_vert /1,2,3,4, 5,6,7,8, 1,2,6,5,
     .                 2,3,7,6, 3,4,8,7, 1,4,8,5/
c
      data nface_axes /1,2, 1,2, 1,3, 2,3, 1,3, 2,3/
c
      data xi /0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0,
     .         1.d0,1.d0,0.d0, 0.d0,1.d0,0.d0,
     .         0.d0,0.d0,1.d0, 1.d0,0.d0,1.d0,
     .         1.d0,1.d0,1.d0, 0.d0,1.d0,1.d0/
c
      data pi /3.14159265358979312d0/
c
      theta_max = pi/2.d0*7.d0/9.d0
c
c
      iprint=0
      allocate(nvisit(NRPOINT))
      nvisit(1:NRPOINT) = 0
c
      nrsplit=0
      do ih=1,NRHEXAS
        do iv=1,8
          np = HEXAS(ih)%VertNo(iv)
          if (POINTS(np)%Type.eq.'CoorNrm') then
            nflag(iv)=1
          else
            nflag(iv)=0
          endif
        enddo
c
c  .....check if there are not more than one face on the reconstructed
c       surface
        ii=0
        do if=1,6
          call decode(HEXAS(ih)%FigNo(if), nr,lab)
          idec=1
          do ifv=1,4
            iv = nface_vert(ifv,if)
            if (nflag(iv).eq.0) idec=0
          enddo
c
c  .......check for neighbors, there should be only one
          if ((RECTANGLES(nr)%BlockNo(1).ne.0).and.
     .        (RECTANGLES(nr)%BlockNo(2).ne.0)) idec=0
          if (idec.eq.1) then
            ii=ii+1
            iface(ii) = if
          endif
        enddo
        if (ii.gt.1) then
          write(*,8001) ih,0,ii
 8001     format(i6,2x,6i3)
          nrsplit=nrsplit+1
        endif
      enddo
      if (nrsplit.gt.0) then
        write(*,8002) nrsplit
 8002   format('check_reconstruct: NUMBER OF HEXAS TO SPLIT = ',i4)
        stop 1
      endif
c
c  ...loop through the hexas
      do ih=1,NRHEXAS
        if (iprint.eq.5) then
          write(*,5038) ih, HEXAS(ih)%VertNo(1:8)
 5038     format('check_reconstruct: ih, VERTICES = ',
     .           i5,3x,8i5)
        endif
c
c  .....determine vertex coordinates
        do iv=1,8
          np = HEXAS(ih)%VertNo(iv)
          xv(1:3,iv) = POINTS(np)%Rdata(1:3)
          if (POINTS(np)%Type.eq.'CoorNrm') then
            nflag(iv)=1
            rn(1:3,iv) = POINTS(np)%Rdata(4:6)
          else
            nflag(iv)=0
            rn(1:3,iv) = 0.d0
          endif
        enddo
c
c  .....check determinants at vertices
        do iv=1,8
          call shapeHb(xi(1:3,iv),norder, nrdof,vshap,dvshap)
          do i=1,3
            dxdxi(1:3,i,iv) = 0.d0
            do k=1,nrdof
              dxdxi(1:3,i,iv) = dxdxi(1:3,i,iv)
     .                        + xv(1:3,k)*dvshap(i,k)
            enddo
          enddo
          call mixed_product(dxdxi(1:3,1,iv),dxdxi(1:3,2,iv),
     .                       dxdxi(1:3,3,iv), rjac)
          if (iprint.eq.1) then
            write(*,7004) ih,iv,rjac
 7004       format('check_reconstruct: ih,iv,rjac = ',
     .                  i5,i2,e12.5)
          endif
          if (rjac.lt.0.d0) then
            write(*,7001) ih,iv
 7001       format('check_reconstruct: NEGATIVE JACOBIAN at ih,iv = ',
     .             i5,i2)
            call pause
          endif
        enddo
c
c  .....check normals for faces
        do if=1,6
          idec=1
          do ifv=1,4
            iv = nface_vert(ifv,if)
            if (nflag(iv).eq.0) idec=0
          enddo
          if (idec.eq.1) then
c
c  .........face on the reconstructed surface
c
c  .........loop through the face vertices
            do ifv=1,4
              iv = nface_vert(ifv,if)
              np = HEXAS(ih)%VertNo(iv)
c
c  ...........compute components of the normal wrt local vertex
c             coordinates
              do k=1,3
                ivk = neig_vert(k,iv)
                aa(1:3,k) = xv(1:3,ivk) - xv(1:3,iv)
              enddo
              bb(1:3) = POINTS(np)%Rdata(4:6)
              if (iprint.eq.5) then
                do i=1,3
                  write(*,7027) i,aa(i,1:3),bb(i)
 7027             format('i,aa,bb = ',i2,2x,3f8.3,2x,f8.3)
                enddo
              endif
              call decomp(3,3,aa,ipiv,iflag)
              if (iflag.eq.1) then
                write(*,*) 'check_reconstruct: SINGULAR MATRIX'
                stop 1
              endif
              call gauss2(3,3,aa,ipiv,bb, uu)
              if (iprint.eq.5) then
                write(*,7007) np,ih,iv,if,ifv,uu
 7007           format('check_reconstruct: np,ih,iv,if,ifv,uu = ',
     .                 i5,i5,3i2,2x,3f8.3)
                call pause
              endif
c
              if ((uu(1).ge.0.d0).and.(uu(2).ge.0.d0).and.
     .          (uu(3).ge.0.d0)) then
                if (nvisit(np).eq.0) then
                  write(*,7005) np,ih,if,ifv
 7005             format('check_reconstruct: MODIFYING NORMAL FOR ',
     .                   'np,ih,if,ifv =',i5,i5,2i2)
                  nvisit(np)=1
                  POINTS(np)%Rdata(4:6) = -POINTS(np)%Rdata(4:6)
                else
                  write(*,7006) np,ih,if,ifv
 7006             format('check_reconstruct: CANNOT REVERSE NORMAL ',
     .                   'FOR np,ih,if,ifv =',i5,i5,2i2)
                  stop 1
                endif
              endif
            enddo
          endif
        enddo
c
c  ...end of loop through hexas
      enddo
      deallocate(nvisit)
ccc      return
c
c-----------------------------------------------------------------------
c
c  ...geometrical corrections...........................................
 10   continue
c
c  ...number of corrected normals
      normals_corrected=0
c
c  ...loop through the hexas
      do ih=1,NRHEXAS
c
c  .....determine vertex coordinates
        do iv=1,8
          np = HEXAS(ih)%VertNo(iv)
          xv(1:3,iv) = POINTS(np)%Rdata(1:3)
          if (POINTS(np)%Type.eq.'CoorNrm') then
            nflag(iv)=1
          else
            nflag(iv)=0
          endif
        enddo
c
c  .....check normals for faces
        do if=1,6
          idec=1
          do ifv=1,4
            iv = nface_vert(ifv,if)
            if (nflag(iv).eq.0) idec=0
          enddo
          if (idec.eq.1) then
c
c  .........face on the reconstructed surface
c
c  .........loop through the face vertices
            do ifv=1,4
              iv = nface_vert(ifv,if)
              np = HEXAS(ih)%VertNo(iv)
c
c  ...........enforce minimum angle condition
c
c  ...........look for the vertex connected to the face vertex
c             that it is NOT on the face
              do k=1,3
                ivk = neig_vert(k,iv)
                call locate(ivk,nface_vert(1:4,if),4, ii)
                if (ii.eq.0) exit
              enddo
c
c  ...........edge vector
              temp(1:3) = xv(1:3,iv) - xv(1:3,ivk)
              call normalize(temp)
              call scalar_product(temp,POINTS(np)%Rdata(4:6), s)
              if (s.lt.cos(theta_max)) then
c
c  .............modify the component along the edge vector
                temp1(1:3) = POINTS(np)%Rdata(4:6) - s*temp(1:3)
                call normalize(temp1)
                POINTS(np)%Rdata(4:6) = sin(theta_max*.99d0)*temp1(1:3)
     .                                + cos(theta_max*.99d0)*temp(1:3)
                write(*,7047) np,ih,POINTS(np)%Rdata(1:3)
 7047           format('input_reconstruct: CORRECTED NORMAL AT np,ih =',
     .                  1x,2i5,3f8.3)
                normals_corrected = normals_corrected+1
              endif
            enddo
          endif
        enddo
c
c  ...end of loop through hexas
      enddo
      if (normals_corrected.gt.0) then
        write(*,*) 'normals_corrected = ',normals_corrected
        call pause
        go to 10
      endif
c
c
      end


