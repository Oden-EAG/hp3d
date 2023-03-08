c----------------------------------------------------------------------
c
c   routine name       - input_RECONSTRUCT
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 11
c
c   purpose            - routine reads in geometry in the RECONSTRUCT
c                        format:
c
c.....................................................................
c
c     NRSURFS:  number of surfaces to be reconstructed
c               'RecSurf' repeated NRSURFS times
c     NRDOMAIN: number of subdomains
c     NRPOINTS: number of points
c               points coordinates
c     NRTETRA:  number of tetrahedra
c               for each tetrahedron:
c                 subdomain number, vertices
c               number of G^1 (or surface) triangles
c               for each G^1 triangle:
c                 surface number, vertices
c
c.....................................................................
c
c                        The routine generates then a full GMP
c                        information, defining curves and remaining
c                        triangles and constructing curve_to_point
c                        and triangle_to_point connectivities.
c                        The full GMP connectivities are completed
c                        in 'connect'. The triangles on reconstructed
c                        surfaces undergo the geometry reconstruction.
c
c---------------------------------------------------------------------
c
      subroutine input_RECONSTRUCT(Fp)
c
      use control
      use GMP
      use element_data
      use pml
c
#include "syscom.blk"
c
      integer, parameter          :: nin = 15
      character(len=*),intent(in) :: Fp
c
c  ...point to tets connectivities (temporary)
      integer, dimension(:,:), pointer :: point_to_tets
c
c  ...number of tets connected to a point (temporary)
      integer, dimension(:),   pointer :: point_nrtets
c
      integer, dimension(2,20) :: list
      integer, dimension(3)    :: nvert_aux
c     for i-th block on the list
c       list(1,i) = block connected to an edge
c       list(2,i) = the corresponding edge number
c     or
c       list(1,i) = block connected to a face
c       list(2,i) = the corresponding face number
c
c     nvert_aux - vertices for a surface triangle
      character(10) surfstype(10), surfrtype(10)
c
c  ...list of surfaces containing a curve
      integer lists(10)
c
c  ...workspace
      dimension rdata_save(3),rpos(3)
c---------------------------------------------------------------------
c
      iprint=0
c
c  ...check user defined parameters NDIM and MANDIM
      if (NDIM.ne.3) then
        write(*,1000)NDIM
1000    format(' input_RECONSTRUCT: inconsistent NDIM = ',i1)
        NDIM = 3
        write(*,1001)NDIM
1001    format('   resetting NDIM = ',i1)
      endif
      if (MANDIM.ne.3) then
        write(*,1002)MANDIM
1002    format(' input_RECONSTRUCT: inconsistent MANDIM = ',i1)
        MANDIM = 3
        write(*,1003)MANDIM
1003    format('   resetting MANDIM = ',i1)
      endif
c
c  ...set up the unused parameters
      NRRECTA = 0;  NRHEXAS = 0;  NRPRISM = 0;  NRPYRAM = 0
c
c  ...allocate memory for GMP data structure
      call alloc_GMP
      write(*,*) 'GMP: data structure allocated'
c
c---------------------------------------------------------------------
c
c  Step 1: Input data on reconstructed surfaces
c
      open(unit=nin,file=Fp,form='formatted',
     .     access='sequential',status='unknown')
c
c  ...read in number of surfaces
      read(nin,*) NRSURFS
      if (MAXSU.lt.NRSURFS) then
        write(*,7010)
 7010   format('input_RECONSTRUCT: INCREASE MAXSU')
        stop
      endif
      if (iprint.eq.1) then
        write(*,7011)
 7011   format('input_RECONSTRUCT: READING SURFACES...')
      endif
c
      do ns=1,NRSURFS
        read(nin,*) SURFACES(ns)%Type
        if (iprint.eq.1) then
          write(*,5001) ns,SURFACES(ns)%Type
 5001     format('input_RECONSTRUCT: ns = ',i3,' SURFACE TYPE = ',a10)
        endif
c
        select case(SURFACES(ns)%Type)
c
c  .......a plane normal to a given vector and passing through
c         a point, orientation specified by the normal
          case('VecPt')
            allocate(SURFACES(ns)%Rdata(6))
            read(nin,*) SURFACES(ns)%Rdata(1:3)
            read(nin,*) SURFACES(ns)%Rdata(4:6)
c
c  .......plane passing through three given points, A,B,C
c         orientation specified by AB x AC
          case('ThrPt')
            allocate(SURFACES(ns)%Rdata(9))
            read(nin,*) SURFACES(ns)%Rdata(1:3)
            read(nin,*) SURFACES(ns)%Rdata(4:6)
            read(nin,*) SURFACES(ns)%Rdata(7:9)
c
c  .......a sphere with radius and central point
          case('Sphere')
            allocate(SURFACES(ns)%Rdata(4))
            read(nin,*) SURFACES(ns)%Rdata(1:3)
            read(nin,*) SURFACES(ns)%Rdata(4)
c
c    .....a cylinder defined with an endpoint of its axis, direction
c         vector and radius
          case('Cylinder')
            allocate(SURFACES(ns)%Rdata(7))
            read(nin,*) SURFACES(ns)%Rdata(1:3)
            read(nin,*) SURFACES(ns)%Rdata(4:6)
            read(nin,*) SURFACES(ns)%Rdata(7)
c
c  .......ellipsoid with semi-axes a,b,c along x,y,z and given central
c         point
          case('Ellipsoid')
            allocate(SURFACES(ns)%Rdata(6))
c
c  .........read coordinates of the central point
            read(nin,*) SURFACES(ns)%Rdata(1:3)
c
c  .........read semiaxes lengths
            read(nin,*) SURFACES(ns)%Rdata(4:6)
c
c  .......cone
          case('Cone')
            allocate(SURFACES(ns)%Rdata(7))
c
c  .........read coordinates of the central point
            read(nin,*) SURFACES(ns)%Rdata(1:3)
c
c  .........read direction vector
            read(nin,*) SURFACES(ns)%Rdata(4:6)
c
c  .........read the cone angle
            read(nin,*) SURFACES(ns)%Rdata(7)
c
c  .......reconstructed surface
          case('RecSurf')
c
c  .......other cases
          case default
            write(*,7012)
 7012       format('input_RECONSTRUCT: UNKNOWN SURFACE TYPE')
            write(*,*) SURFACES(ns)%Type
            stop 1
        end select
      enddo
c
c  ...input number of subdomains
      read(nin,*) NRDOMAIN
      write(*,*) 'GMP: Surface and domain read'
c
c---------------------------------------------------------------------
c
c  Step 2: Input points
c
c
c  ...read in number of points
      read(nin,*) NRPOINT
      if (MAXNP.lt.NRPOINT) then
        write(*,*) 'input_RECONSTRUCT:INCREASE MAXNP'
        stop 1
      endif
      if (iprint.eq.1) then
        write(*,*) 'input_RECONSTRUCT: READING POINTS...'
        write(*,*) '              NRPOINT = ',NRPOINT
      endif
      do np=1,NRPOINT
        POINTS(np)%Type  = 'Regular'
        allocate(POINTS(np)%Rdata(3))
        read(nin,*) POINTS(np)%Rdata(1:3)
      enddo
      write(*,*) 'GMP: Points read'
c
c---------------------------------------------------------------------
c
c  Step 3: Input tetrahedra
c
c  ...read in number of tets
      read(nin,*) NRTETRA
      if (NRTETRA.gt.MAXTE) then
        write(*,7013) NRTETRA, MAXTE
 7013   format('input_RECONSTRUCT: INCREASE MAXTE,NRTETRA,MAXTE = ',2i8)
        stop 1
      endif
      do ntet=1,NRTETRA
        TETRAS(ntet)%Type = 'Linear'
        TETRAS(ntet)%FigNo(1:4) = 0
        TETRAS(ntet)%EdgeNo(1:6) = 0
        read(nin,*) TETRAS(ntet)%domain,TETRAS(ntet)%VertNo(1:4)
c
c  .....check orientation
        ! call check_orientation(3,ntet)
      enddo
      write(*,*) 'GMP: Tetras read'
c
c---------------------------------------------------------------------
c
c  Step 4: Determine point to tets connectivities
c
c  ...max anticipated number of tetras connected to a point
      mptet=140
      allocate(point_to_tets(1:mptet,1:NRPOINT))
      allocate(point_nrtets(1:NRPOINT))
      point_nrtets(1:NRPOINT) = 0
c
c  ...loop through tets
      do ntet=1,NRTETRA
c
c  .....loop through tetrahedron vertex points
        do jv=1,4
c
          np = TETRAS(ntet)%VertNo(jv)
          nrtets = point_nrtets(np)
c
c  .......add the tetrahedron to the list of tetras connected to the
c         point
          nrtets = nrtets+1
          if (nrtets.gt.mptet) then
            write(*,7021) mptet
 7021       format('input_RECONSTRUCT: increase mptet = ',i3)
            stop 1
          endif
          point_to_tets(nrtets,np) = ntet
          point_nrtets(np) = nrtets
c
c  .....end loop through vertices
        enddo
c
c  ...end of loop through tetras
      enddo
      if (iprint.eq.1) then
        do np=1,NRPOINT
          write(*,7022) np, (point_to_tets(ii,np),ii=1,point_nrtets(np))
 7022     format('input_RECONSTRUCT: np = ',i4,' CONNECTED TETS = ',
     .           /,10i8,/,10i8,/,10i8,/,10i8,/,10i8,/,10i8,/,10i8/,10i8)
        enddo
        call pause
      endif
      write(*,*) 'GMP: Point to Tet connectivity'
c
c-----------------------------------------------------------------------
c
c  Step 5: Count curves and determine temporary tets to curves
c          connectivities
c
      NRCURVE=0
c
c  ...loop through tetras
      do ntet=1,NRTETRA
c
c  .....loop through the tetrahedron edge curves
        do je=1,6
c
c  .......if the edge curve has not been defined and connected yet
          if (TETRAS(ntet)%EdgeNo(je).eq.0) then
c
c  .........initiate the counter for tets that are
c           connected to the edge
            ibl=0
c
c  .........determine the endpoints of the edge
            iv1 = TETRA_EDGE_TO_VERT(1,je)
            np1 = TETRAS(ntet)%VertNo(iv1)
            iv2 = TETRA_EDGE_TO_VERT(2,je)
            np2 = TETRAS(ntet)%VertNo(iv2)
c
c  .........loop through the tetras connected to the first endpoint
            nrtets = point_nrtets(np1)
            do k=1,nrtets
              ntet1 = point_to_tets(k,np1)
c
c  ...........locate the vertex of 'ntet1' coinciding with the first
c             endpoint
              call locate(np1,TETRAS(ntet1)%VertNo(1:4),4, kv1)
c
c  ...........lower the flag indicating that 'ntet1' shares the edge
              iflag=0
c
c  ...........loop through the edges of 'ntet1' sharing the vertex
              do l=1,3
                ie = TETRA_VERT_TO_EDGE(l,kv1)
c
c  .............look for the other vertex point
                do ive=1,2
                  iv3 = TETRA_EDGE_TO_VERT(ive,ie)
                  np3 = TETRAS(ntet1)%VertNo(iv3)
c
c  ...............raise the flag if the other endpoint has been found
                  if (np2.eq.np3) iflag=1
                enddo
c
c  .............stop the loop if the edge has been found
                if (iflag.eq.1) exit
c
c  ...........end of loop through the edges sharing the vertex
              enddo
c
c  ...........'ntet1' shares the point but not the edge, skip it
              if (iflag.eq.0) go to 10
c
c  ...........pick up the edge curve number for ntet1
              nc = TETRAS(ntet1)%EdgeNo(ie)
c
c  ...........if the curve has already been counted
              if (nc.ne.0) then
                write(*,7023) ntet,je,ntet1
 7023       format('input_RECONSTRUCT: INCONSISTENCY ntet,je,ntet1 = ',
     .             i8,i3,i8)
                stop
              else
c
c  .............add 'tet1' to the list of tets to be
c               connected to the curve
                ibl=ibl+1
                list(1,ibl) = ntet1; list(2,ibl) = ie
              endif
c
 10           continue
c
c  .........end of loop through connected blocks
            enddo
c
c  .........create a new curve
            NRCURVE = NRCURVE+1
            if (NRCURVE.gt.MAXNC) then
              write(*,7024) MAXNC
 7024         format('input_RECONSTRUCT: INCREASE MAXNC = ',i10)
              stop 1
            endif
            nc = NRCURVE
            if (iprint.eq.1) then
              write(*,7025) nc,ntet,je
 7025         format('input_RECONSTRUCT: HAVE DEFINED CURVE ',i8,
     .               ' FOR TET ', i8,' AND EDGE ',i2)
              write(*,*) 'CONNECTED TETS WITH EDGE NUMBERS'
              write(*,8001) (list(1,jj),jj=1,ibl)
 8001         format(10i8)
              write(*,8001) (list(2,jj),jj=1,ibl)
ccc              call pause
            endif
c
c  .........connect all blocks from the list to the newly created curve
            nrbl = ibl
            do ibl=1,nrbl
              nbl = list(1,ibl)
              ie  = list(2,ibl)
              TETRAS(nbl)%EdgeNo(ie) = nc
            enddo
c
c  .......if the edge has not been connected yet
          endif
c
c  .....end of loop through edges of the tet
        enddo
c
c  ...end of loop through tetras
      enddo
c
c  ...define curves..................
      do nc=1,NRCURVE
        CURVES(nc)%Type = 'void'
      enddo
c
c  ...loop through tets
      do ntet=1,NRTETRA
c
c  .....loop through the tetra edge curves
        do je=1,6
c
c  .......get the connected curve
          nc = TETRAS(ntet)%EdgeNo(je)
c
c  .......if the curve is visited for the first type
          if (CURVES(nc)%Type.eq.'void') then
c
c  .........determine the endpoints of the edge
            iv1 = TETRA_EDGE_TO_VERT(1,je)
            np1 = TETRAS(ntet)%VertNo(iv1)
            iv2 = TETRA_EDGE_TO_VERT(2,je)
            np2 = TETRAS(ntet)%VertNo(iv2)
            CURVES(nc)%Type = 'Seglin'
            CURVES(nc)%EndPoNo(1) = np1
            CURVES(nc)%EndPoNo(2) = np2
            CURVES(nc)%NrFig = 0
          endif
c
c  .....end of loop through edges of the tetrahedron
        enddo
c
c  ...end of loop through tetras
      enddo
c
c  ...erase the temporary block to edge curves connectivities,
c     they will be regenarated in routine connect accounting
c     for orientation
      do ntet=1,NRTETRA
        TETRAS(ntet)%EdgeNo(1:6) = 0
      enddo
      if (iprint.eq.1) then
        do nc=1,NRCURVE
          write(*,7026) nc,(CURVES(nc)%EndPoNo(ii),ii=1,2)
 7026     format('input_RECONSTRUCT: nc,points = ',i8,2x,2i8)
        enddo
        call pause
      endif
      write(*,*) 'GMP: Curve to Tet connectivity'
c
c----------------------------------------------------------------------
c
c  Step 6: count triangles and determine temporary tets to face
c          triangles connectivities
c
      NRTRIAN=0
c
c  ...loop through tets
      do ntet=1,NRTETRA
c
c  .....loop through the tetra face triangles
        do jf=1,4
c
c  .......if the face has not been connected yet
          if (TETRAS(ntet)%FigNo(jf).eq.0) then
c
c  .........initiate counter for the tetras connected to the face
            ibl=0
c
c  .........determine vertex points on the face
            iv = TETRA_FACE_TO_VERT(1,jf)
            np1 = TETRAS(ntet)%VertNo(iv)
            iv = TETRA_FACE_TO_VERT(2,jf)
            np2 = TETRAS(ntet)%VertNo(iv)
            iv = TETRA_FACE_TO_VERT(3,jf)
            np3 = TETRAS(ntet)%VertNo(iv)
c
c  .........loop through the tetras connected to the first vertex point
            nrbl = point_nrtets(np1)
            do k=1,nrbl
              ntet1 = point_to_tets(k,np1)
c
c  ...........locate the vertex coinciding with the point
              call locate(np1,TETRAS(ntet1)%VertNo(1:4),4, iv1)
c
c  ...........loop through the faces sharing the point
              do l=1,3
                if = TETRA_VERT_TO_FACE(l,iv1)
c
c  .............look for two other common points
                iflag=0
                do ivf=1,3
                  iv = TETRA_FACE_TO_VERT(ivf,if)
                  np = TETRAS(ntet1)%VertNo(iv)
                  if (np.eq.np2) iflag=iflag+1
                  if (np.eq.np3) iflag=iflag+1
                enddo
                if (iflag.eq.2) exit
              enddo
              if (iflag.ne.2) go to 30
c
c  ...........look for a connected triangle
              nt = TETRAS(ntet1)%FigNo(if)
c
c  ...........if the triangle has already been created
              if (nt.ne.0) then
                write(*,7031) ntet,if,ntet1
 7031        format('input_RECONSTRUCT: INCONSISTENCY ntet,if,ntet1 = ',
     .              i8,i3,i8)
                stop
              else
c
c  .............add the block to the list of blocks
c               to be connected to the curve
                ibl=ibl+1
                list(1,ibl) = ntet1; list(2,ibl) = if
              endif
c
 30           continue
c
c  .........end of loop through connected blocks
            enddo
c
c  .........no block connected to the face has been found, the
c           face rectangle has to be created
            NRTRIAN = NRTRIAN + 1
            if (NRTRIAN.gt.MAXTR) then
              write(*,7032) MAXTR
 7032         format('input_RECONSTRUCT: INCREASE MAXTR = ',i10)
              stop 1
            endif
            nt = NRTRIAN
            if (iprint.eq.1) then
              write(*,7033) nt,ntet,jf
 7033         format('input_RECONSTRUCT: HAVE DEFINED TRIANGLE ',i8,
     .               ' FOR TETRA ', i8,' AND FACE ',i2)
              write(*,*) 'CONNECTED TETS WITH FACE NUMBERS'
              write(*,8001) (list(1,jj),jj=1,ibl)
              write(*,8001) (list(2,jj),jj=1,ibl)
ccc              call pause
            endif
c
c  .........connect all blocks from the list to the rectangle
            nrbl = ibl
            do ibl=1,nrbl
              nn = list(1,ibl)
              if = list(2,ibl)
              TETRAS(nn)%FigNo(if) = nt
            enddo
c
c  .......if the face has not been connected yet
          endif
c
c  .....end of loop through faces of the tetra
        enddo
c
c  ...end of loop through tetras
      enddo
      if (iprint.eq.1) call pause
c
c
c  ...define the triangles
      do nr=1,NRTRIAN
        TRIANGLES(nr)%Type = 'void'
      enddo
c
c  ...loop through tetra
      do ntet=1,NRTETRA
c
c  .....loop through the tetra face triangles
        do jf=1,4
          nt = TETRAS(ntet)%FigNo(jf)
c
c  .......if visited for the first time
          if (TRIANGLES(nt)%Type.eq.'void') then
c
c  .........determine vertex points on the face
            iv = TETRA_FACE_TO_VERT(1,jf)
            np1 = TETRAS(ntet)%VertNo(iv)
            iv = TETRA_FACE_TO_VERT(2,jf)
            np2 = TETRAS(ntet)%VertNo(iv)
            iv = TETRA_FACE_TO_VERT(3,jf)
            np3 = TETRAS(ntet)%VertNo(iv)
            TRIANGLES(nt)%Type = 'PlaneTri'
            TRIANGLES(nt)%VertNo(1) = np1
            TRIANGLES(nt)%VertNo(2) = np2
            TRIANGLES(nt)%VertNo(3) = np3
            TRIANGLES(nt)%BlockNo(1:2) = 0
ccc            TRIANGLES(nt)%Surface = 0
c
          endif
c
c  .....end of loop through faces of the tetra
        enddo
c
c  ...end of loop through tetras
      enddo
c
c  ...erase the temporary block to face triangles connectivities,
c     they will be reconnected in routine connect accounting
c     for orientation
      do ntet=1,NRTETRA
        TETRAS(ntet)%FigNo(1:4) = 0
      enddo
      if (iprint.eq.1) then
        do nt = 1,NRTRIAN
          write(*,7034) nt,(TRIANGLES(nt)%VertNo(ii),ii=1,3)
 7034     format('input_GMPdata: nt,points = ',i8,2x,3i8)
        enddo
ccc        call pause
      endif
      write(*,*) 'GMP: Trian to Tet'
c
c----------------------------------------------------------------------
c
c  Step 7: for each algebraic or reconstructed surface, read in the surface
c          triangles, recheck the compatibility of connected tetrahedra
c          and redefine the triangle
c
c  ...read the number of surface triangles
      read(nin,*) no_surface_triangles
c
c  ...loop through surface triangles
      do itsurf=1,no_surface_triangles
c
c  .....read in data
        read(nin,*) nos, nvert_aux(1:3)
c
        if (IHEAD_FLAG.eq.0) then
c
c     .....find the corresponding triangle in the GMP data structure
          call find_GMPtriangle(nvert_aux(1:3), nt)
          if (nt.eq.0) then
            write(*,7041) nos, nvert_aux(1:3)
 7041       format('input_NETGEN: HAVE NOT FOUND A TRIANGLE FOR ',
     .        'SURFACE ',i2,' VERTICES = ',3i6)
            call print_GMP
            stop
          endif
c
c     .....store surface number in data structure
          select case (SURFACES(nos)%Type)
c     .....reconstructed surface triangle
          case('RecSurf')
ccc   TRIANGLES(nt)%Type = 'G1RecTri'
            TRIANGLES(nt)%Type = 'PlaneTri'
            allocate(TRIANGLES(nt)%Idata(1))
            TRIANGLES(nt)%Idata(1) = nos
c     .....algebraic surface, use parametric transfinite interpolation
          case default
            TRIANGLES(nt)%Type = 'PTITri'
            allocate(TRIANGLES(nt)%Idata(1))
            TRIANGLES(nt)%Idata(1) = nos
          end select
        endif
c
c  ...end of loop through surface triangles
      enddo
c
      if (IHEAD_FLAG.eq.1) then
        write(*,*) 'GMP: IHEAD_FLAG on: manual setting for head problem'
        do nr=1,NRTRIAN
          rpos = POINTS(TRIANGLES(nr)%VertNo(1))%Rdata(1:3)
          r1   = sqrt(rpos(1)**2 + rpos(2)**2 + rpos(3)**2)
          rpos = POINTS(TRIANGLES(nr)%VertNo(2))%Rdata(1:3)
          r2   = sqrt(rpos(1)**2 + rpos(2)**2 + rpos(3)**2)
          rpos = POINTS(TRIANGLES(nr)%VertNo(3))%Rdata(1:3)
          r3   = sqrt(rpos(1)**2 + rpos(2)**2 + rpos(3)**2)
c
          rsum = abs(r1-RPML_MIN)+abs(r2-RPML_MIN)+abs(r3-RPML_MIN)
          if ((rsum.lt.GEOM_TOL)) then
            TRIANGLES(nr)%Type = 'PTITri'
            allocate(TRIANGLES(nr)%Idata(1))
            TRIANGLES(nr)%Idata(1) = PML_SPHERE
          endif
        enddo
      endif
c
      if (iprint.eq.10) then
        write(*,7046)
 7046   format('input_RECONSTRUCT: HAVE REDEFINED SURFACE TRIANGLES')
      endif
c  ...deallocate temporary connectivities
      deallocate(point_to_tets,point_nrtets)
c
      write(*,*) 'GMP: Temporary workspace deallocated'
c
c
!=======================================================================
! G E O M E T R Y    I N F O    S O    F A R :                         |
!                                                                      |
!        POINTS%Type         = 'Regular'                               |
!              %Rdata(1:3)   = coordinates                             |
!                                                                      |
!        CURVES%Type         = 'Seglin'                                |
!              %EndPoNo(1:2) = endpoints numbers                       |
!                                                                      |
!        TRIANGLES%Type      = 'PlaneTri','PTITri'                     |
!                 %VertNo    = vertex numbers                          |
!                 %Idata(1)  = surface number, if 'PTITri' or          |
!                              'PlaneTri' on a 'RecSurf'               |
!                                                                      |
!        TETRAS%Type         = 'Linear'                                |
!              %VertNo(1:4)  = vertex numbers                          |
!              %Domain       = domain number                           |
!                                                                      |
!=======================================================================
c
      call upgrade2double
      write(*,*) 'GMP: upgrade2double finished'
c
!=======================================================================
!  IF NEEDED, GEOMETRY CUSTOMIZATION MUST BE DONE HERE                 |
!=======================================================================
c
      call customize_geometry
      write(*,*) 'GMP: customization'
c-----------------------------------------------------------------------
c
c  STEP 9: complete GMP connectivities
c
      call connect
      write(*,*) 'GMP: connected'
c
c-----------------------------------------------------------------------
c
c  STEP 10: redefine curvilinear segments on ALGEBRAIC surfaces
c
      do nc = 1, NRCURVE
c
c  .....initiate the number of different algebraic and reconstructed
c       surface triangles adjacent to the curve segment
        i = 0; its = 0; itr = 0
c
c  .....loop through figures connected to the curve
        do if=1,CURVES(nc)%NrFig
          nick = abs(CURVES(nc)%FigNo(if))
          call decode(nick, nf,lab)
          select case(lab)
c
c  .........triangle
            case(1)
ccc            if ((TRIANGLES(nf)%Type.eq.'G1RecTri').or.
ccc     .          (TRIANGLES(nf)%Type.eq.'PTITri')) then
            if (TRIANGLES(nf)%Type.eq.'PTITri') then
c
c  ...........find the surface number
              ns = TRIANGLES(nf)%Idata(1)
              call locate(ns, lists,i, num)
              if (num.eq.0) then
                i=i+1; lists(i) = ns
                if (SURFACES(ns)%Type.eq.'RecSurf') then
                  itr=itr+1
                  surfrtype(itr) = 'G1RecTri'
                else
                  its=its+1
                  surfstype(its) = SURFACES(ns)%Type
                endif
              endif
            endif
c
c  .........rectangle
            case(2)
ccc            if ((RECTANGLES(nf)%Type.eq.'G1RecRec').or.
            if (RECTANGLES(nf)%Type.eq.'PTIRec') then
c
c  ...........find the surface number
              ns = RECTANGLES(nf)%Idata(1)
              call locate(ns, lists,i, num)
              if (num.eq.0) then
                i=i+1; lists(i) = ns
                if (SURFACES(ns)%Type.eq.'RecSurf') then
                  itr=itr+1
                  surfrtype(itr) = 'G1RecRec'
                else
                  its=its+1
                  surfstype(its) = SURFACES(ns)%Type
                endif
              endif
            endif
           end select
        enddo
c
c  .....curve on a single algebraic surface
        if ((its.eq.1).and.(itr.eq.0)) then
          CURVES(nc)%Type = '1SurfsCur'
          allocate(CURVES(nc)%Idata(1))
          CURVES(nc)%Idata(1) = lists(1)
c
c  .....intersection of two algebraic surfaces
        elseif ((its.eq.2).and.(itr.eq.0)) then
          CURVES(nc)%Type = '2SurfsCur'
          allocate(CURVES(nc)%Idata(2))
          CURVES(nc)%Idata(1:2) = lists(1:2)
c
c======================================================================
c  REMARK: do nothing for reconstructed surfaces.                     |
c                                                                     |
c  .....curve on a single reconstructed surface                       |
        elseif ((its.eq.0).and.(itr.eq.1)) then
ccc          CURVES(nc)%Type = '1SurfrCur'
ccc          allocate(CURVES(nc)%Idata(1))
ccc          CURVES(nc)%Idata(1) = lists(1)
cccc
c  .....intersection of two reconstructed surfaces (sharp edge)
        elseif ((its.eq.0).and.(itr.eq.2)) then
ccc          CURVES(nc)%Type = '2SurfrCur'
ccc          allocate(CURVES(nc)%Idata(2))
ccc          CURVES(nc)%Idata(1:2) = lists(1:2)
c======================================================================
c
c  .....intersection of three algebraic surfaces
        elseif ((its .eq. 3) .and. (itr .eq. 0)) then
          CURVES(nc)%Type = '3SurfsCur'
          allocate(CURVES(nc)%Idata(3))
          CURVES(nc)%Idata(1:3) = lists(1:3)
c
c  .....intersection of four algebraic surfaces
        elseif ((its .eq. 4) .and. (itr .eq. 0)) then
          CURVES(nc)%Type = '4SurfsCur'
          allocate(CURVES(nc)%Idata(4))
          CURVES(nc)%Idata(1:4) = lists(1:4)
c
c  .....other cases not supported for now...
        elseif ((its.ne.0).or.(itr.ne.0)) then
          write(*,7051) nc,its,itr
 7051     format('input_RECONSTRUCT: nc,its,itr = ',i7,2i3,
     .           ' CURVE CASE NOT SUPPORTED')
          write(*,7052) surfstype(1:its)
 7052     format('              SURFACE TYPES: ',10(2x,a10))
          call print_GMP
        endif
c
c  ...end of loop through curves
      enddo
c
c----------------------------------------------------------------------
c
c  STEP 11: redefine curvilinear triangles and rectangles
c
      do nt=1,NRTRIAN
        if (TRIANGLES(nt)%Type.eq.'PlaneTri') then
          nflag=0
          do ie=1,3
            nc = abs(TRIANGLES(nt)%EdgeNo(ie))
            if (CURVES(nc)%Type.ne.'Seglin') nflag=1
          enddo
          if (nflag.eq.1) TRIANGLES(nt)%Type = 'TransTri'
        endif
      enddo
c
      do nr=1,NRRECTA
        if (RECTANGLES(nr)%Type.eq.'BilQua') then
          nflag=0
          do ie=1,4
            nc = abs(RECTANGLES(nr)%EdgeNo(ie))
            if (CURVES(nc)%Type.ne.'Seglin') nflag=1
          enddo
          if (nflag.eq.1) RECTANGLES(nr)%Type = 'TraQua'
        endif
      enddo
c
c----------------------------------------------------------------------
c
c  STEP 12: redefine curvilinear tets, prisms and pyramids
c
      do ntet=1,NRTETRA
        nflag=0
        do ie=1,6
          nc = abs(TETRAS(ntet)%EdgeNo(ie))
          if (CURVES(nc)%Type.ne.'Seglin') nflag=1
        enddo
        do if=1,4
          if (TETRAS(ntet)%FigNo(if).eq.0) cycle
          call decode(TETRAS(ntet)%FigNo(if), nt,lab)
          if (TRIANGLES(nt)%Type.ne.'PlaneTri') nflag=1
        enddo
        if (nflag.eq.1) then
          if (iprint.eq.1) then
            write(*,*) 'input_RECONSTRUCT: REDEFINING ntet = ',ntet
          endif
          TETRAS(ntet)%Type = 'TraTet'
        endif
      enddo
c
      do npri=1,NRPRISM
        nflag=0
        do ie=1,9
          nc = abs(PRISMS(npri)%EdgeNo(ie))
          if (CURVES(nc)%Type.ne.'Seglin') nflag=1
        enddo
        do if=1,2
          if (PRISMS(npri)%FigNo(if).eq.0) cycle
          call decode(PRISMS(npri)%FigNo(if), nt,lab)
          if (TRIANGLES(nt)%Type.ne.'PlaneTri') nflag=1
        enddo
        do if=3,5
          if (PRISMS(npri)%FigNo(if).eq.0) cycle
          call decode(PRISMS(npri)%FigNo(if), nr,lab)
          if (RECTANGLES(nr)%Type.ne.'BilQua') nflag=1
        enddo
        if (nflag.eq.1) then
          if (iprint.eq.1) then
            write(*,*) 'input_RECONSTRUCT: REDEFINING npri = ',npri
          endif
          PRISMS(npri)%Type = 'TIprism'
        endif
      enddo
c
      do npyr=1,NRPYRAM
        nflag=0
        do ie=1,8
          nc = abs(PYRAMIDS(npyr)%EdgeNo(ie))
          if (CURVES(nc)%Type.ne.'Seglin') nflag=1
        enddo
        do if=1,1
          if (PYRAMIDS(npyr)%FigNo(if).eq.0) cycle
          call decode(PYRAMIDS(npyr)%FigNo(if), nr,lab)
          if (RECTANGLES(nr)%Type.ne.'BilQua') nflag=1
        enddo
        do if=2,5
          if (PYRAMIDS(npyr)%FigNo(if).eq.0) cycle
          call decode(PYRAMIDS(npyr)%FigNo(if), nt,lab)
          if (TRIANGLES(nt)%Type.ne.'PlaneTri') nflag=1
        enddo
        if (nflag.eq.1) then
          if (iprint.eq.1) then
            write(*,*) 'input_RECONSTRUCT: REDEFINING npyr = ',npyr
          endif
          PYRAMIDS(npyr)%Type = 'TIpyram'
        endif
      enddo
c
c======================================================================
c  SO FAR, RECONSTRUCTED SURFACES HAVE BEEN DEFINED AS PIECEWISE      |
c  LINEAR                                                             |
c======================================================================
c
c----------------------------------------------------------------------
c
      close(nin)
c
      if (iprint.eq.1) write(*,*)'input_RECONSTRUCT: done!'
c
      write(*,*) 'GMP: Reconstruction done'
c
      return










c
c  STEP 13: redefine points on reconstructed surfaces
c
c  ...loop over 'G1RecTri'
      do i = 1,NRTRIAN
        if (TRIANGLES(i)%Type.ne.'G1RecTri') cycle
        do iv = 1,3
          np = TRIANGLES(i)%VertNo(iv)
c  .......if vertex has already been visited, cycle
          if (POINTS(np)%Type.eq.'CoorNrm') cycle
          POINTS(np)%Type = 'CoorNrm'
          rdata_save(1:3) = POINTS(np)%Rdata(1:3)
          deallocate(POINTS(np)%Rdata, STAT=is)
          if (is.ne.0) then
       write(*,*)'input_RECONSTRUCT: Rdata not deallocated for np = ',np
            stop
          endif
          allocate(POINTS(np)%Rdata(6), STAT=is)
          if (is.ne.0) then
       write(*,*)'input_RECONSTRUCT: Rdata not allocated for np = ',np
            stop
          endif
          POINTS(np)%Rdata(1:3) = rdata_save(1:3)
        enddo
c  ...end of loop over 'G1RecTri'
      enddo

c----------------------------------------------------------------------
c
c  STEP 14: redefine curves type and reconstruct G^1 geometry
c
      do i=1,NRCURVE
        select case(CURVES(i)%Type)
        case('1SurfrCur','2SurfrCur')
          CURVES(i)%Type='HermCur'
        endselect
      enddo
c
      close(nin)
c
      if (iprint.eq.1) write(*,*)'input_RECONSTRUCT: done!'
c
c-----------------------------------------------------------------------
c
      call reconstruct
ccc 100  continue
ccc      write(*,*)'Enter domain to be smoothed (0 - exit)'
ccc      read(*,*)ndom
ccc      if (ndom.gt.0) then
ccc        call smooth_domain(ndom)
ccc        goto 100
ccc      endif
c
c
      end subroutine input_RECONSTRUCT
