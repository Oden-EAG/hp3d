c----------------------------------------------------------------------
c
c   routine name       - input_trianG1
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 10
c
c----------------------------------------------------------------------
c
      subroutine input_trianG1
c
      use control
      use GMP
      use element_data
c      
#include "syscom.blk"     
#include "cinout.blk"
c
      integer, dimension(:),   pointer :: vpoint_to_spoint
      integer, dimension(:,:), pointer :: point_to_blocks
      integer, dimension(:),   pointer :: point_nrbl
c
c  ...work space
      dimension rdata(6)
c
c  ...maximum anticipated number of tets attached to a point      
      parameter (mphex = 20)
c
c  ...work space      
      integer list(2,20)
c     for i-th block on the list
c       list(1,i) = block connected to an edge
c       list(2,i) = the corresponding edge number
c     or
c       list(1,i) = block connected to a face
c       list(2,i) = the corresponding face number 
c
      character*10 type
c
c----------------------------------------------------------------------
c
      iprint=0
c      
      write(*,*)'entered input_trianG1...'
c
c  ...set up the dimension of the problem
      NDIM=3
c
c  ...set up the unused parameters
      NRSURFS=0; NRRECTA=0; NRPRISM=0; NRHEXAS=0; NRPYRAM=0
c
c----------------------------------------------------------------------
c
c  Step 1: Input points
c
c======================================================================      
c  REMARK: the numbering of the points is dictated by the order of    |
c          the SURFACE points, ie surface points are read first, then |
c          volume points are read and, if volume point is not on the  |
c          surface, a new point is created.                           |
c======================================================================      
c
c  ...read the number of points on the surface and in the
c     TOTAL (!) number of points, hence including surface points
      read(KIN,*) nrsp, NRPOINT
c
c  ...allocate the memory dynamically
      allocate(POINTS(NRPOINT), STAT=is) 
      if (is.ne.0) then
        write(*,*) 'input_trianG1: POINTS ARE NOT ALLOCATED'
        stop
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'input_trianG1: READING SURFACE POINTS...'
      endif
c
c**********************************************************************
c  ...read in the surface points
      do np=1,nrsp
        POINTS(np)%Type  = 'CoorNrm'
        POINTS(np)%NrCurv = 0
        allocate(POINTS(np)%Rdata(6), STAT=is)
        if (is.ne.0) then 
          write(*,*) 'input_trianG1: Rdata IS  NOT ALLOCATED',
     .               'POINT NUMBER = ', np
          stop
        endif
        read(KIN,*) POINTS(np)%Rdata(1:6)
c  .....renormalize the normal
        call normalize(POINTS(np)%Rdata(4:6))
      enddo   
c
c  ...allocate volume points to surface points connectivities
      allocate(vpoint_to_spoint(NRPOINT), STAT=is) 
      if (is.ne.0) then
        write(*,*) 'input_trianG1: vpoint_to_spoint ALLOCATION ',
     .             ' FAILED'
        stop
      endif
      vpoint_to_spoint(1:NRPOINT)=0
c
      if (iprint.eq.1) then
        write(*,*) 'input_trianG1: READING VOLUME POINTS'
      endif
c
c**********************************************************************
c  ...read in the volume points
      np=nrsp
c  ...loop through volume points      
      do ip=1,NRPOINT
        read(KIN,*) rdata(1:3)
c  .....search for a corresponding point on the surface
        do jp=1,nrsp
          dist=0.d0
          do ivar=1,3
            dist = dist + (POINTS(jp)%Rdata(ivar)-rdata(ivar))**2
          enddo
          dist = sqrt(dist)
c  .......have found a matching point on the surface
          if (dist.lt.GEOM_TOL) then
c  .........save the connectivity
            vpoint_to_spoint(ip)=jp
          endif
c  .....end of loop through surface points
        enddo
c
c  .....if no surface point has been found, ie point is INSIDE 
c       the surface
        if (vpoint_to_spoint(ip).eq.0) then
c  .......generate a new regular point
          np=np+1
          if (np.gt.NRPOINT) then
            write(*,8001) np,NRPOINT
 8001       format('input_trianG1: INCONSISTENCY 1, np,NRPOINT = ',
     .             2i6)
            stop
          endif
          POINTS(np)%Type  = 'Regular'
          POINTS(np)%NrCurv = 0
          allocate(POINTS(np)%Rdata(3), STAT=is)
          if (is.ne.0) then 
            write(*,*) 'input_trianG1: Rdata IS  NOT ALLOCATED',
     .                 'POINT NUMBER = ', np
            stop
          endif
          POINTS(np)%Rdata(1:3) = rdata(1:3)
c  .......save the connectivity
          vpoint_to_spoint(ip)=np
        endif
c  ...end of loop through volume points
      enddo

      if (iprint.ge.1) then
        write(*,*) 'input_trianG1: HAVE READ POINTS'
      endif
c
c----------------------------------------------------------------------
c
c  Step 2: Input tets
c
c  ...read the number of hexas in the volume
      read(KIN,*) NRTETRA
      if (iprint.eq.1) then
        write(*,*) 'input_trianG1: NRTETRA = ',NRTETRA
      endif
c
c  ...set up the dimension of the manifold
      if (NRTETRA.le.0) then
        write(*,8002) NRTETRA
 8002   format('input_trianG1: NRTETRA = ',i6)
        stop
      else
        MANDIM=3
      endif
c
c  ...allocate dynamically the memory for the tets
      allocate(TETRAS(NRTETRA), STAT=is ) 
      if (is.ne.0) then
        write(*,*) 'input_trianG1: TETS ARE NOT ALLOCATED'
        stop
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'input_trianG1: READING TETS...'
      endif
c 
c  ...loop over tets      
      do ih=1,NRTETRA
c  .....set tet's type equal to transfinite interpolation      
        TETRAS(ih)%Type = 'TraTet'
        read(KIN,*) TETRAS(ih)%VertNo(1:4)
c  .....change the enumeration if needed
        do i=1,4
          j = TETRAS(ih)%VertNo(i)
          TETRAS(ih)%VertNo(i) = vpoint_to_spoint(j)
        enddo
        if (iprint.eq.1) then
          write(*,7636) ih, TETRAS(ih)%VertNo(1:4)
 7636     format('HEXA = ',i6,' VERTICES = ',8i7)
        endif
c  ...end of loop through tets        
      enddo
c
c  ...printing statement      
      if (iprint.ge.1) then
        write(*,*) 'input_trianG1: HAVE READ TETS'
        call pause
      endif
c
c----------------------------------------------------------------------
c
c  Step 3: Determine point to tets connectivities
c
c  ...allocate work spaces      
      allocate(point_to_blocks(1:mphex,1:NRPOINT), STAT=is ) 
      allocate(point_nrbl(1:NRPOINT), STAT=is)
      point_nrbl(1:NRPOINT) = 0
c
c  ...loop through tets
      do nh=1,NRTETRA    
c  .....compute nickname for tet      
        nick = nh*10+3
c  .....loop through tet vertex points
        do jv=1,4
          np = TETRAS(nh)%VertNo(jv)
          nrbl = point_nrbl(np)
c  .......check if the tet is on the list of blocks connected to the point
          call locate(nick,point_to_blocks(1:nrbl,np),nrbl, ii)
c  .......if tet was not found          
          if (ii.eq.0) then
c  .........update number of connected tets                  
            nrbl = nrbl+1
            if (nrbl.gt.mphex) then
              write(*,*) 'input_trianG1: INCREASE mphex = ',mphex
              stop
            endif
c  .........store tet's nickname            
            point_to_blocks(nrbl,np) = nick
            point_nrbl(np) = nrbl 
          endif
c  .....end loop through vertices
        enddo
c  ...end of loop through tets
      enddo

c  ...printing statements      
      if (iprint.eq.1) then
        do np=1,NRPOINT
          write(*,7005) np, (point_to_blocks(ii,np),ii=1,point_nrbl(np))
 7005     format('input_trianG1: np = ',i7,
     .           ' CONNECTED BLOCKS = ',15i8)
        enddo
        call pause
      endif
      if (iprint.ge.1) then
        write(*,*) 'input_trianG1: HAVE CONNECTED TETS TO POINTS'
      endif
c
c----------------------------------------------------------------------
c
c  Step 4a: Count curves and determine temporary connectivities
c
      NRCURVE=0
c
c  ...initialize edge number to 0      
      do nh=1,NRTETRA
        TETRAS(nh)%EdgeNo(1:6) = 0
      enddo
c
c
c
c PG, Feb 10: NEEDS TO BE FINISHED!!!!       
c     LOOP over tets 
c     | LOOP over edges
c     | | IF 1st visit to edge
c     | | | get 1st and 2nd endpoint
c     | | | LOOP over tets connected to 1st endpoint
c     | | | | locate local number of endpoint
c     | | | | LOOP over edges connected to endpoint
c     | | | | | search for 2nd endpoint
c     | | | | | IF found EXIT
c     | | | | ENDLOOP
c     | | | | IF not found GOTO 10
c     | | | | get curve number
c     | | | | IF not 1st visit
c     | | | | | store curve number
c     | | | | | GOTO 20
c     | | | | ELSE
c     | | | | | add tet to list
c     | | | | | GOTO 10
c     | | | | ENDIF
c  10
c     | | | ENDLOOP
c     | | | add curve
c
cc
c
c  ...loop through tets
      do nh=1,NRTETRA  
c  .....compute nickname for tet      
        nick = nh*10+3
c  .....loop through the tet edge curves
        do je=1,6  
c  .......if the edge has not been connected yet
          if (TETRAS(nh)%EdgeNo(je).eq.0) then
c  .........initiate the list of blocks connected to the edge
c           with the current block
            ibl=1
            list(1,ibl) = nick
            list(2,ibl) = je
c  .........determine the endpoints of the edge
            iv1 = TETRA_EDGE_TO_VERT(1,je)
            np1 = TETRAS(nh)%VertNo(iv1)
            iv2 = TETRA_EDGE_TO_VERT(2,je)
            np2 = TETRAS(nh)%VertNo(iv2)
c  .........loop through the tets connected to 1st endpoint
            nrbl = point_nrbl(np1)
            do k=1,nrbl
              nick1 = point_to_blocks(k,np1)
              call decode(nick1, nh1,lab)
c  ...........select block type              
              select case(lab)
c  .............tetrahedron
                case(3)
c  ...............locate the vertex coinciding with the point
                  call locate(np1,TETRAS(nh1)%VertNo(1:4),4, kv1)
                  iflag=0
c  ...............loop through the edges sharing the vertex
                  do l=1,3
                    ie = TETRA_VERT_TO_EDGE(l,kv1)
c  .................look for the other vertex point
                    do ive=1,2
                      iv3 = TETRA_EDGE_TO_VERT(ive,ie)
                      np3 = TETRAS(nh1)%VertNo(iv3)
                      if (np2.eq.np3) iflag=1
                    enddo
                    if (iflag.eq.1) exit
c  ...............end of loop through the edge sharing the vertex
                  enddo
                  if (iflag.eq.0) go to 10
                  nc = TETRAS(nh1)%EdgeNo(ie)
c  ...............if the curve has already been counted
                  if (nc.ne.0) then
c  .................store the connection for the original hexa
                    TETRAS(nh)%EdgeNo(je) = nc  
                    go to 20
                  else
c  .................add the block to the list of blocks
c                   to be connected to the curve
                    ibl=ibl+1
                    list(1,ibl) = nick1; list(2,ibl) = ie
c  .................proceed to the next connected block
                    go to 10
                  endif
                case default
                  write(*,*) 'input_trianG1: WRONG lab = ',lab
                  stop
c  ...........end select block type                  
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
 7003         format('input_trianG1: HAVE DEFINED CURVE ',i4,
     .               ' FOR TET ', i4,' AND EDGE ',i2)
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
c  ...........tetrahedron              
              case(3)
                ie = list(2,ibl)
                TETRAS(nn)%EdgeNo(ie) = nc
              end select
            enddo
c
c  .......if the edge has not been connected yet
          endif
 20       continue
c
c  .....end of loop through edges of the tet
        enddo
c
c  ...end of loop through tets
      enddo
c
c  ...define curves..................
      allocate(CURVES(NRCURVE), STAT=is) 
      if (is.ne.0) then
        write(*,*) 'input_trianG1: CURVES ARE NOT ALLOCATED'
        stop
      endif
      do nc=1,NRCURVE
        CURVES(nc)%Type = 'void'
      enddo
c
c  ...loop through tets
      do nh=1,NRTETRA  
c
c  .....loop through the tet's edge curves
        do je=1,6  
c
c  .......get the connected curve
          nc = TETRAS(nh)%EdgeNo(je)
c
c  .......if the curve is visited for the first type
          if (CURVES(nc)%Type.eq.'void') then
            if (iprint.eq.1) then
              write(*,*) 'input_trianG1: DEFINING CURVE nc = ',nc
            endif
c
c  .........determine the endpoints of the edge
            iv1 = TETRA_EDGE_TO_VERT(1,je)
            np1 = TETRAS(nh)%VertNo(iv1)
            iv2 = TETRA_EDGE_TO_VERT(2,je)
            np2 = TETRAS(nh)%VertNo(iv2)
            CURVES(nc)%Type = 'Seglin'
            CURVES(nc)%EndPoNo(1) = np1
            CURVES(nc)%EndPoNo(2) = np2
            CURVES(nc)%NrFig = 0
          endif
c
c  .....end of loop through edges of the tet
        enddo
c
c  ...end of loop through tets
      enddo
c
c  ...erase the temporary block to edge curves connectivities
      do nh=1,NRTETRA
        TETRAS(nh)%EdgeNo(1:6) = 0
      enddo
      if (iprint.eq.1) then
        do nc = 1,NRCURVE
          write(*,7001) nc,(CURVES(nc)%EndPoNo(ii),ii=1,2)
 7001     format('input_trianG1: nc,points = ',i4,2x,2i5)
        enddo
        call pause
      endif
      if (iprint.ge.1) then
        write(*,*) 'input_trianG1: HAVE DEFINED AND CONNECTED',
     .             ' CURVES'
      endif
c
c----------------------------------------------------------------------
c
c  Step 5: count triangles and determine temporary tets to face triangle
c          connectivities
c
      NRTRIAN=0
c
c  ...initialize fig number to 0      
      do nh=1,NRTETRA
        TETRAS(nh)%FigNo(1:6) = 0
      enddo
c
c  ...loop through tetra
      do nh=1,NRTETRA  
c
c  .....store nickname for tet      
        nick = nh*10+3
c
c  .....loop through the tet faces
        do jf=1,4
c
c  .......if the face has not been connected yet
          if (TETRAS(nh)%FigNo(jf).eq.0) then
c
c  .........initiate the list of tets connected to the face
c           to the current block
            ibl=0
c
c  .........determine global number of face vertices
            iv = TETRA_FACE_TO_VERT(1,jf)
            np1 = TETRAS(nh)%VertNo(iv)
            iv = TETRA_FACE_TO_VERT(2,jf)
            np2 = TETRAS(nh)%VertNo(iv)
            iv = TETRA_FACE_TO_VERT(3,jf)
            np3 = TETRAS(nh)%VertNo(iv)
c
c  .........determine number of tets connected to the first vertex
            nrbl = point_nrbl(np1)
c  .........loop over tets connected to the first vertex
            do k=1,nrbl
c  ...........determine tet's nickname            
              nick1 = point_to_blocks(k,np1)
c  ...........determine tet's number              
              call decode(nick1, nh1,lab)
c              
c  ...........select element type              
              select case(lab)
c
c  ...........tetrahedron
              case(3)
c  .............determine local number of vertex, when visited
c               from tet nh1
                call locate(np1,TETRAS(nh1)%VertNo(1:4),4, iv1)
c
c  .............loop through the faces sharing the point 
                do l=1,3
c  ...............determine local number of face                
                  if = TETRA_VERT_TO_FACE(l,iv1)
c
c  ...............look for two other common points
                  iflag=0
                  do ivf=1,3
                    iv = TETRA_FACE_TO_VERT(ivf,if)
                    np = TETRAS(nh1)%VertNo(iv)
                    if (np.eq.np2) iflag=iflag+1
                    if (np.eq.np3) iflag=iflag+1
                  enddo
                  if (iflag.eq.2) exit
                enddo
                if (iflag.ne.2) go to 30
c
c  .............look for a connected triangle
                nf = TETRAS(nh1)%FigNo(if)
c
c  .............if the triangle has already been created
                if (nf.ne.0) then
c
c  ...............store the connection for the original tet
                  TETRAS(nh)%FigNo(jf) = nf
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
                write(*,*) 'input_trianG1: WRONG lab = ',lab
                stop
              end select
c
 30           continue
c
c  .........end of loop through connected blocks
            enddo
c
c  .........no block connected to the face has been found, the
c           face rectangle has to be created
            NRTRIAN = NRTRIAN + 1
            nr = NRTRIAN
c  .........printing statement            
            if (iprint.eq.5) then
              write(*,7008) nr,nh,jf
 7008         format('input_trianG1: HAVE DEFINED RECTA ',i4,
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
c  ...........select block type              
              select case(lab)
c  ...........tetrahedron              
              case(3)
                if = list(2,ibl)
                TETRAS(nn)%FigNo(if) = nr
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
      allocate( TRIANGLES(NRTRIAN), STAT=is ) 
      if (is.ne.0) then
        write(*,*) 'input_trianG1: TRIANGLES ARE NOT ALLOCATED'
        stop
      endif
      do nr=1,NRTRIAN
        TRIANGLES(nr)%Type = 'void'
      enddo
c
c  ...loop through tets
      do nh=1,NRTETRA  
c
c  .....loop through the tet's faces
        do jf=1,4  
          nr = TETRAS(nh)%FigNo(jf)
c
c  .......if visited for the first time
          if (TRIANGLES(nr)%Type.eq.'void') then
c  .........determine vertex points on the face
            iv = TETRA_FACE_TO_VERT(1,jf)
            np1 = TETRAS(nh)%VertNo(iv)
            iv = TETRA_FACE_TO_VERT(2,jf)
            np2 = TETRAS(nh)%VertNo(iv)
            iv = TETRA_FACE_TO_VERT(3,jf)
            np3 = TETRAS(nh)%VertNo(iv)
            TRIANGLES(nr)%Type = 'PlaneTri'
            TRIANGLES(nr)%VertNo(1) = np1
            TRIANGLES(nr)%VertNo(2) = np2
            TRIANGLES(nr)%VertNo(3) = np3
            TRIANGLES(nr)%BlockNo(1:2) = 0
          endif
c
c  .....end of loop through faces of the hexa
        enddo
c
c  ...end of loop through tets
      enddo
c
c  ...erase the temporary block to face rectangles connectivities
      do nh=1,NRTETRA
        TETRAS(nh)%FigNo(1:4) = 0
      enddo
      if (iprint.eq.1) then
        do nr = 1,NRTRIAN
          write(*,7002) nr,(TRIANGLES(nr)%VertNo(ii),ii=1,3)
 7002     format('input_trianG1: nr,points = ',i4,2x,3i5)
        enddo
        call pause
      endif
c
c      
c----------------------------------------------------------------------
c      
c  Step 6: set up the domain number
c
      NRDOMAIN = 0
      do i=1,NRTETRA
        NRDOMAIN = NRDOMAIN + 1
        TETRAS(i)%Domain = NRDOMAIN
      enddo
c
c
c----------------------------------------------------------------------
c      
c  Step 7: complete connectivities
c
      deallocate(point_to_blocks,point_nrbl, STAT=is)
      deallocate(vpoint_to_spoint, STAT=is)
      if (iprint.ge.1) then
        write(*,*) 'input_trianG1: calling connect'
      endif
      call connect
c
c  ...data consistency check
      if (iprint.ge.1) then
        write(*,*) 'input_trianG1: CHECKING..'
      endif
ccc      call check_reconstruct
c
c
c----------------------------------------------------------------------
c
c  Step 8: update curve type for curves on the surface
c
c  ...loop over curves
      do nc = 1,NRCURVE
        iflag = 1
c  .....loop over curve's endpoints        
        do iv = 1,2
          nv = CURVES(nc)%EndPoNo(iv)
c  .......lower flag if endpoint is not on surface          
          if (POINTS(nv)%Type.ne.'CoorNrm') iflag = 0
c  .....end loop over endpoints          
        enddo
c  .....update curve's type if needed        
        if (iflag.eq.1) CURVES(nc)%Type='cBezier'
c  ...end loop over curves        
      enddo
c      
c----------------------------------------------------------------------
c
c  Step 9: update triangle type for triangles on the surface
c
c  ...loop over triangles
      do nt = 1,NRTRIAN
        iflag = 0
c  .....loop over triangle's edge  
        do ie = 1,3
          nc = abs(TRIANGLES(nt)%EdgeNo(ie))
          if (CURVES(nc)%Type.eq.'cBezier') iflag = iflag + 1
        enddo
c  .....update triangle's type if needed        
        if (iflag.ge.1) TRIANGLES(nt)%Type = 'TransTri'
        if (iflag.eq.3) TRIANGLES(nt)%Type = 'G1RecTri'
c  ...end loop over triangles
      enddo




ccc  ...loop over triangles
cc      do nt = 1,NRTRIAN
cc        iflag = 1
ccc  .....loop over triangle's vertices     
cc        do iv = 1,3
cc          nv = TRIANGLES(nt)%VertNo(iv)
ccc  .......lower flag if endpoint is not on surface          
cc          if (POINTS(nv)%Type.ne.'CoorNrm') iflag = 0
ccc  .....end loop over vertices          
cc        enddo
ccc  .....update triangle's type if needed        
cc        if (iflag.eq.1) TRIANGLES(nt)%Type='G1RecTri'
ccc  ...end loop over triangles
cc      enddo
       
c----------------------------------------------------------------------
c      
c  Step 9: reconstruct the geometry
c
      if (iprint.ge.1) then
ccc        write(*,*) 'input_trianG1: BEGINNING RECONSTRUCTION'
      endif
ccc      call reconstruct
      if (iprint.ge.1) then
ccc        write(*,*) 'input_trianG1: FINISHED RECONSTRUCTION'
      endif
c
c
      write(*,*)'exited input_trianG1...'
c      
      end
