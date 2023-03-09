c----------------------------------------------------------------------
c
c   routine name       - input_reconstrucT
c
c----------------------------------------------------------------------
c
c   latest revision    - May 2023
c
c   purpose: routine performs G1 reconstruction of a triangulated
c            surface
c
c   STRUCTURE
c      input_reconstruct
c      | POINTS: 'Regular','CoorNor'
c      | CURVES: 'Seglin'
c      | HEX
c      |--> reconstruct
c           | RECTANGLES: 'HermRec'
c           | CURVES:     'HermCur'
c           |--> set_tangents
c           |    |--> recontruct_normal
c           |--> reconstruct_curve
c           |--> set_tangents
c           |--> reconstruct_patch
c                mixed derivatives!
c                |--> reconstruct_curvature
c      recta_HermRec


c----------------------------------------------------------------------
c
      subroutine input_reconstruct
c
      use control
      use GMP
      use element_data
c
      implicit none
c
c  ...temporary connectivities workspaces
      integer, dimension(:),   pointer :: vpoint_to_spoint
      integer, dimension(:,:), pointer :: point_to_blocks
      integer, dimension(:),   pointer :: point_nrbl
c
c  ...tet workspace
      integer, dimension(10) :: temp
      real(8), dimension(3)  :: rdata_save
c
c  ...maximum anticipated number of tets attached to a point
      integer, parameter :: mphex = 100
c
c  ...work space
      integer :: list(2,20)
c     for i-th block on the list
c       list(1,i) = block connected to an edge
c       list(2,i) = the corresponding edge number
c     or
c       list(1,i) = block connected to a face
c       list(2,i) = the corresponding face number
c
      character(10) :: type
c
      integer :: iprint
c
c----------------------------------------------------------------------
c
      iprint=0
c
c  ...set up the dimension of the problem
      NDIM=3
c  ...set up the unused parameters
      NRSURFS=0; NRRECTA=0; NRPRISM=0; NRHEXAS=0; NRPYRAM=0
      call alloc_GMP
c
c----------------------------------------------------------------------
c
c  Step 1: input points
c
c  ...read in number of points
      read(KIN,*) NRPOINT
      if (iprint.eq.1) then
        write(*,*)'input_reconstruct: NRPOINT = ',NRPOINT
      endif
      if (NRPOINT.gt.MAXNP) then
        write(*,1000)
 1000   format(' input_reconstruct: INCREASE MAXNP!')
        stop
      endif

c  ...allocate memory for points
ccc      allocate(POINTS(NRPOINT), STAT=is)
ccc      if (is.ne.0) then
ccc        write(*,*)'input_reconstruct: POINTS not allocated!'
ccc        stop
ccc      endif
c  ...printing
      if (iprint.eq.1) then
        write(*,*)'input_reconstruct: reading points...'
      endif
c  ...loop over points
      do np = 1,NRPOINT
        POINTS(np)%Type   = 'Regular'
        POINTS(np)%NrCurv = 0
        allocate(POINTS(np)%Rdata(3), STAT=is)
        if (is.ne.0) then
          write(*,*) 'input_reconstruct: Rdata not allocated!',
     .               'np = ', np
          stop
        endif
        read(KIN,*) POINTS(np)%Rdata(1:3)
        if (iprint.eq.1) then
          write(*,*)'i,Rdata = ',np,POINTS(np)%Rdata(1:3)
        endif
c  ...end of loop over points
      enddo
c  ...printing
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: have read points.'
        call pause
      endif
c
c----------------------------------------------------------------------
c
c  Step 2: input tets
c
c  ...read the number of tetras in the volume
      read(KIN,*) NRTETRA
      if (iprint.eq.1) then
        write(*,*)'input_reconstruct: NRTETRA = ',NRTETRA
      endif
      if (NRTETRA.gt.MAXTE) then
        write(*,1001)
 1001   format(' input_reconstruct: INCREASE MAXTE!')
        stop
      endif

c  ...read the number of domains
      read(KIN,*) NRDOMAIN
      if (iprint.eq.1) then
        write(*,*)'input_reconstruct: NRDOMAIN = ',NRDOMAIN
      endif
c  ...set up the dimension of the manifold
      if (NRTETRA.le.0) then
        write(*,8002) NRTETRA
 8002   format('input_reconstruct: NRTETRA = ',i6)
        stop
      else
        MANDIM=3
      endif
c  ...allocate memory for the tets
ccc      allocate(TETRAS(NRTETRA), STAT=is)
ccc      if (is.ne.0) then
ccc        write(*,*) 'input_reconstruct: TETRAS not allocated!'
ccc        stop
ccc      endif
c  ...printing
      if (iprint.eq.1) then
        write(*,*) 'input_reconstruct: reading tets...'
      endif
c  ...loop over tets
      do n = 1,NRTETRA
c  .....get tet's data
        read(KIN,*)TETRAS(n)%Domain,TETRAS(n)%VertNo(1:4)
        TETRAS(n)%Type = 'Linear'
        if (iprint.eq.1) then
          write(*,7636) n,TETRAS(n)%VertNo(1:4),TETRAS(n)%Domain
 7636     format('TET = ',i6,' VERTICES = ',8i7,' DOMAIN = ',i3)
        endif
c  ...end of loop through tets
      enddo
c  ...check tets orientation and make it positive if needed
      do n = 1,NRTETRA
        if (iprint.eq.1) then
          write(*,*)'checking orientation for n = ',n
        endif
        call check_orientation(3,n)
      enddo
c  ...printing
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: have read tets.'
        call pause
      endif
c
c----------------------------------------------------------------------
c
c  Step 3: determine point to tets connectivities
c
c  ...allocate work spaces
      allocate(point_to_blocks(1:mphex,1:NRPOINT), STAT=is )
      allocate(point_nrbl(1:NRPOINT), STAT=is)
      point_nrbl(1:NRPOINT) = 0
c
c  ...loop through tets
      do nh = 1,NRTETRA
c  .....compute nickname for tet
        nick = nh*10 + 3
c  .....loop through tet vertex points
        do jv = 1,4
          np = TETRAS(nh)%VertNo(jv)
          nrbl = point_nrbl(np)
c  .......check if the tet is on the list of blocks connected to the point
          call locate(nick,point_to_blocks(1:nrbl,np),nrbl, ii)
c  .......if tet was not found
          if (ii.eq.0) then
c  .........update number of connected tets
            nrbl = nrbl + 1
            if (nrbl.gt.mphex) then
              write(*,*) 'input_reconstruct: INCREASE mphex = ',mphex
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
 7005     format('input_reconstruct: np = ',i7,
     .           ' CONNECTED BLOCKS = ',15i8)
        enddo
        call pause
      endif
      if (iprint.ge.1) then
        write(*,*) 'input_reconstruct: have connected tets to point'
      endif
c
c----------------------------------------------------------------------
c
c  Step 4a: count curves and determine temporary connectivities
c
      NRCURVE = 0
c
c  ...initialize edge number to 0
      do nh = 1,NRTETRA
        TETRAS(nh)%EdgeNo(1:6) = 0
      enddo
c
c  ...loop through tets
      do nh = 1,NRTETRA
c  .....compute nickname for tet
        nick = nh*10 + 3
c  .....loop through the tet edge curves
        do je = 1,6
c  .......if the edge has not been connected yet
          if (TETRAS(nh)%EdgeNo(je).eq.0) then
c  .........initiate the list of blocks connected to the edge
c           with the current block
            ibl = 1
            list(1,ibl) = nick
            list(2,ibl) = je
c  .........determine the endpoints of the edge
            iv1 = TETRA_EDGE_TO_VERT(1,je)
            np1 = TETRAS(nh)%VertNo(iv1)
            iv2 = TETRA_EDGE_TO_VERT(2,je)
            np2 = TETRAS(nh)%VertNo(iv2)
c  .........loop through the tets connected to 1st endpoint
            nrbl = point_nrbl(np1)
            do k = 1,nrbl
              nick1 = point_to_blocks(k,np1)
              call decode(nick1, nh1,lab)
c  ...........select block type
              select case(lab)
c  .............tetrahedron
                case(3)
c  ...............locate the vertex coinciding with the point
                  call locate(np1,TETRAS(nh1)%VertNo(1:4),4, kv1)
                  iflag = 0
c  ...............loop through the edges sharing the vertex
                  do l = 1,3
                    ie = TETRA_VERT_TO_EDGE(l,kv1)
c  .................look for the other vertex point
                    do ive = 1,2
                      iv3 = TETRA_EDGE_TO_VERT(ive,ie)
                      np3 = TETRAS(nh1)%VertNo(iv3)
                      if (np2.eq.np3) iflag = 1
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
                  write(*,*) 'input_reconstruct: WRONG lab = ',lab
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
            NRCURVE = NRCURVE + 1
            if (NRCURVE.gt.MAXNC) then
              write(*,*)'input_reconstruct: increase MAXNC ',NRCURVE
              stop
            endif
            nc = NRCURVE
            if (iprint.eq.1) then
              write(*,7003) nc,nh,je
 7003         format('input_reconstruct: HAVE DEFINED CURVE ',i4,
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
            do ibl = 1,nrbl
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
ccc      allocate(CURVES(NRCURVE), STAT=is)
ccc      if (is.ne.0) then
ccc        write(*,*) 'input_reconstruct: CURVES not allocated!'
ccc        stop
ccc      endif
      do nc = 1,NRCURVE
        CURVES(nc)%Type = 'void'
      enddo
c
c  ...loop through tets
      do nh = 1,NRTETRA
c
c  .....loop through the tet's edge curves
        do je = 1,6
c
c  .......get the connected curve
          nc = TETRAS(nh)%EdgeNo(je)
c
c  .......if the curve is visited for the first type
          if (CURVES(nc)%Type.eq.'void') then
            if (iprint.eq.1) then
              write(*,*) 'input_reconstruct: DEFINING CURVE nc = ',nc
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
      do nh = 1,NRTETRA
        TETRAS(nh)%EdgeNo(1:6) = 0
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
c  Step 5: count triangles and determine temporary tets to face triangle
c          connectivities
c
      NRTRIAN = 0
c
c  ...initialize fig number to 0
      do nh = 1,NRTETRA
        TETRAS(nh)%FigNo(1:4) = 0
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
                write(*,*) 'input_reconstruct: WRONG lab = ',lab
                stop
              end select
c
 30           continue
c
c  .........end of loop through connected blocks
            enddo
c
c  .........no block connected to the face has been found, the
c           face triangle has to be created
            NRTRIAN = NRTRIAN + 1
            if (NRTRIAN.gt.MAXTR) then
              write(*,*)'input_reconstruct: increase MAXTR',NRTRIAN
              stop
            endif
            nr = NRTRIAN
c  .........printing statement
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
c  .........connect all blocks from the list to the triangle
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
c  .....end of loop through faces of the tet
        enddo
c
c  ...end of loop through tets
      enddo
c
c
c  ...define the triangles......................
ccc      allocate( TRIANGLES(NRTRIAN), STAT=is)
ccc      if (is.ne.0) then
ccc        write(*,*) 'input_reconstruct: TRIANGLES ARE NOT ALLOCATED'
ccc        stop
ccc      endif
      do nr = 1,NRTRIAN
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
c  ...erase the temporary block to face triangle connectivities
      do nh=1,NRTETRA
        TETRAS(nh)%FigNo(1:4) = 0
      enddo
      if (iprint.eq.1) then
        do nr = 1,NRTRIAN
          write(*,7002) nr,(TRIANGLES(nr)%VertNo(ii),ii=1,3)
 7002     format('input_reconstruct: nr,points = ',i4,2x,3i5)
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

      icount = 0
      do i = 1,NRTRIAN
        if (TRIANGLES(i)%BlockNo(2).eq.0) icount = icount + 1
      enddo
      write(*,*)'icount = ',icount
      call pause

c
c======================================================================
c  Activate to check linear geometry                                  |
      call graphg                                                  !
c  Activate to export mesh                                            |
ccc      call export_mesh_phantom('files/phantom')                    !
ccc      stop                                                         !
c======================================================================
c
c----------------------------------------------------------------------
c
c  STEP 7: update geometry definitions to curvilinear
c
c  ...loop over triangles
      do nt = 1,NRTRIAN
c  .....if only 1 adjacent block
        if (TRIANGLES(nt)%BlockNo(2).eq.0) then
c  .......update triangle to 'G1RecTri'
          TRIANGLES(nt)%Type = 'G1RecTri'
c  .......update adjacent tet to 'TraTet'
          nick = TRIANGLES(nt)%BlockNo(1)
          call decode(nick, ntet,lab)
          select case(lab)
c  .......tetrahedron
          case(3)
            TETRAS(ntet)%Type = 'TraTet'
          case default
            write(*,*)'input_reconstruct: case not supported yet!'
            stop
          endselect
c  .......loop over triangle vertices and edges
          do i = 1,3
            nv = TRIANGLES(nt)%VertNo(i)
            ne = abs(TRIANGLES(nt)%EdgeNo(i))
c  .........update vertices to 'CoorNrm'
            POINTS(nv)%Type = 'CoorNrm'
            rdata_save(1:3) = POINTS(nv)%Rdata(1:3)
            deallocate(POINTS(nv)%Rdata, STAT=is)
            if (is.ne.0) then
              write(*,*)'input_reconstruct: Rdata not deallocated!'
              write(*,*)'                   nv = ',nv
              stop
            endif
            allocate(POINTS(nv)%Rdata(6), STAT=is)
            if (is.ne.0) then
              write(*,*)'input_reconstruct: Rdata not allocated!'
              write(*,*)'                   nv = ',nv
              stop
            endif
            POINTS(nv)%Rdata = 0.d0
            POINTS(nv)%Rdata(1:3) = rdata_save(1:3)
c  .........update edge curves to 'HermCur'
            CURVES(ne)%Type = 'HermCur'
          enddo
        endif
c  ...end of loop over triangles
      enddo
c
c  ...loop over triangles
      do nt = 1,NRTRIAN
c  .....if not a 'PlaneTri' triangle cycle
        if (TRIANGLES(nt)%Type.ne.'PlaneTri') cycle
c  .....loop over edges
        do i = 1,3
          nc = abs(TRIANGLES(nt)%EdgeNo(i))
c  .......if edge is 'HermCur', update triangle to 'TransTri'
          if (CURVES(nc)%Type.eq.'HermCur') then
            TRIANGLES(nt)%Type = 'TransTri'
            exit
          endif
c  .....end of loop over edges
        enddo
c  ...end of loop over triangles
      enddo
c
c  ...loop over tetra
      do n = 1,NRTETRA
c  .....loop over tetra edges
        do i = 1,6
          nc = abs(TETRAS(n)%EdgeNo(i))
c  .......update tet type if needed
          if (CURVES(nc)%Type.eq.'HermCur')  TETRAS(n)%Type = 'TraTet'
c  .....end of loop over edges
        enddo
c  ...end of loop over tetra
      enddo
c
c
c======================================================================
c  Activate to export to STL file                                     |
ccc      call export_mesh2STL('files/mesh.stl')                    !
ccc      write(*,*)'mesh exported to STL format!'
ccc      stop
c======================================================================
c
c----------------------------------------------------------------------
c
c  Step 8: reconstruct the geometry
c
      if (iprint.ge.1) then
        write(*,*)'input_reconstruct: BEGINNING RECONSTRUCTION...'
      endif
c
      call reconstruct
      call blow_up(2)
ccc      call thicken_G1_surf_new(2,0.0005d0,0.0005d0)
ccc      call check_jacobian_prism
      if (iprint.ge.1) then
        write(*,*)'input_reconstruct: RECONSTRUCTION SUCCESSFUL!'
      endif
c
      end subroutine input_reconstruct
