!----------------------------------------------------------------------
!> @brief Data structure arrays for Geometry Modeling Package (GMP)
!
!> @date Mar 2023
!----------------------------------------------------------------------
!
   module GMP
!
      implicit none
!
!----------------------------------------------------------------------
!  PARAMETERS                                                         |
!----------------------------------------------------------------------
!     dimension of the problem (2, 3)
      integer, save :: NDIM
!     dimension of the manifold (2, 3)
      integer, save :: MANDIM
!     number of algebraic surfaces
      integer, save :: NRSURFS
!     number of points
      integer, save :: NRPOINT
!     number of curves
      integer, save :: NRCURVE
!     number of triangles
      integer, save :: NRTRIAN
!     number of rectangles
      integer, save :: NRRECTA
!     number of hexahedra
      integer, save :: NRHEXAS
!     number of tetrahedra
      integer, save :: NRTETRA
!     number of prisms
      integer, save :: NRPRISM
!     number of pyramids
      integer, save :: NRPYRAM
!     number of subdomains
      integer, save :: NRDOMAIN
!     number of fine grid points
      integer, save :: NRPO_FG
!     fine grid scale
      real(8), save :: SCALE_FG
!
!     maximum number of algebraic surfaces
      integer, save :: MAXSU
!     maximum number of points
      integer, save :: MAXNP
!     maximum number of curves
      integer, save :: MAXNC
!     maximum number of triangles
      integer, save :: MAXTR
!     maximum number of rectangles
      integer, save :: MAXRE
!     maximum number of prismatic blocks
      integer, save :: MAXBT
!     maximum number of hexahedral blocks
      integer, save :: MAXHE
!     maximum number of tetrahedral blocks
      integer, save :: MAXTE
!     maximum number of pyramids
      integer, save :: MAXPY
!     maximum number of fine grid points
      integer, save :: MAXNP_FG
!
!
!----------------------------------------------------------------------
!  SURFACE                                                            |
!----------------------------------------------------------------------
!
!  Type =
!
!     VecPt     : plane defined by a point and a normal
!       Rdata(1:3) - coordinates of the point
!       Rdata(4:6) - normal
!
!     ThrPt     : plane defined by the coordinates of three points
!       Rdata(1:3) - coordinates of the 1st point
!       Rdata(4:6) - coordinates of the 2nd point
!       Rdata(7:9) - coordinates of the 3rd point
!
!     PtNoVec   : plane defined by a point number and a normal
!       Idata(1)   - point number
!       Rdata(1:3) - normal
!
!     PtNo2Pt   : plane defined by 1 point number and the coordinates
!                 of 2 points
!       Idata(1)   - point number
!       Rdata(1:3) - coordinates of the 1st point
!       Rdata(4:6) - coordinates of the 2nd point
!
!     PPwCC : yz plane parametrized with cylindrical coordinates
!       Rdata(1)   - x coordinate
!       Rdata(2)   - rmin
!       Rdata(3)   - rmax
!
!     Sphere    : sphere
!       Rdata(1:3) - coordinates of the center
!       Rdata(4)   - radius
!
!     Cylinder  : infinite cylinder
!       Rdata(1:3) - coordinates of a point on the axis
!       Rdata(4:6) - axis
!       Rdata(7)   - radius
!
!     Ellipsoid : ellipsoid
!       Rdata(1:3) - coordinates of the center
!       Rdata(4:6) - length of 1st,2nd,3rd semiaxis
!
!     Cone      : cone
!       Rdata(1:3) - coordinates of the vertex
!       Rdata(4:6) - axis
!       Rdata(7)   - half-aperture [Radiants]
!
!     RecSurf   : reconstructed surface
!
!  Idata  - integer data associated to the surface
!
!  Rdata  - real    data associated to the surface
!
!----------------------------------------------------------------------
      type surface
        character(len=10)                       :: Type
        integer,          dimension(:), pointer :: Idata
        real(8),          dimension(:), pointer :: Rdata
      endtype surface
!
!
!----------------------------------------------------------------------
!  POINT                                                              |
!----------------------------------------------------------------------
!
!  Type =
!
!     Regular  : a point
!       Rdata(1:3) - coordinates
!
!     Implicit : implicitly defined point
!       Idata(1:3) - surface nubers identifing the point
!       Rdata(1:3) - coordinates determined by NR iterations
!
!     CoorNrm  : point with a normal
!       Rdata(1:3) - coordinates
!       Rdata(1:3) - normal
!
!     SharpPt  : point on a sharp edge
!       Rdata(1:3) - coordinates
!       Rdata(3:6) - normal to the 1st surface defining the sharp edge
!       Rdata(7:9) - normal to the 2nd surface defining the sharp edge
!
!  NrCurv - number of curves connected to the point
!
!  CurvNo - list of connected curves
!
!  Idata  - integer data associated to the point
!
!  Rdata  - real    data associated to the point
!
!----------------------------------------------------------------------
      type point
        character(len=10)                       :: Type
        integer                                 :: NrCurv
        integer,          dimension(:), pointer :: CurvNo
        integer,          dimension(:), pointer :: Idata
        real(8),          dimension(:), pointer :: Rdata
      endtype point
!
!
!----------------------------------------------------------------------
!  CURVE                                                              |
!----------------------------------------------------------------------
!
!  Type =
!
!     Seglin    : straight segment
!
!     QuaCir    : quarter of a circle
!       Rdata(1:3) - coordinates of the center
!
!     QuaEl1    : quarter of a ellipse (0 to pi/2)
!       Rdata(1:3) - coordinates of the center
!
!     QuaEl2    : quarter of a ellipse (-pi/4 to pi/4)
!       Rdata(1:3) - coordinates of the center
!
!     QuaSEl    : quarter of a superellipse (0 to pi/2)
!       Rdata(1:3) - coordinates of the center
!       Rdata(4:5) - px,py powers of superellipse
!
!     ImpCir    : implicit curve
!       Idata(1:4) - surfaces defining the curve
!
!     3HermCur  : cubic Hermite curve
!       Rdata(1:3) - derivative at 1st end point
!       Rdata(4:6) - derivative at 2nd end point
!
!     1SurfsCur : curve lying on 1 algebraic surface
!       Idata(1) - surface number
!
!     2SurfsCur : curve on the intersection of 2 algebraic surfaces
!       Idata(1:2) - surface numbers
!
!     3SurfsCur : curve on the intersection of 3 algebraic surfaces
!       Idata(1:3) - surface numbers
!
!     HermCur   : quintic Hermite curve
!       Rdata( 1: 3) - derivative at 1st end point
!       Rdata( 4: 6) - derivative at 2nd end point
!       Rdata( 7: 9) - second derivative at 1st end point
!       Rdata(10:11) - second derivative at 2nd end point
!
!     5Bezier   : quintic Bezier curve
!       Rdata(0:17) - control points (6 x 3) from 1st to 2nd endpoint
!
!     7Bezier   : septic Bezier curve
!       Rdata(0:23) - control points (8 x 3) from 1st to 2nd endpoint
!
!     CylCur    : image of a straight line segment through a global
!                 system of cylindrical coordinates
!
!  EndPoNo - curve endpoints numbers
!
!  NrFig   - number of figures connected to the curve
!
!  FigNo   - list of nicknames of connected figures
!     triangle  : nick = ntrian*10 + 1
!     rectangle : nick = nrecta*10 + 2
!
!  Idata  - integer data associated to the curve
!
!  Rdata  - real    data associated to the curve
!
!----------------------------------------------------------------------
      type ccurve
        character(len=10)                       :: Type
        integer                                 :: EndPoNo(2)
        integer                                 :: NrFig
        integer,          dimension(:), pointer :: FigNo
        integer,          dimension(:), pointer :: Idata
        real(8),          dimension(:), pointer :: Rdata
      endtype ccurve
!
!
!----------------------------------------------------------------------
!  TRIANGLE                                                           |
!----------------------------------------------------------------------
!
!  Type =
!
!     PlaneTri : plane triangle
!
!     TransTri : transfinite interpolation triangle
!
!     PTITri   : parametric transfinite interpolation triangle
!       Idata(  1) - conforming surface number
!
!     ImpliTri : implicit triangle
!       Idata(  1) - conforming surface number
!       Idata(2:4) - bounding surfaces numbers (listed counter
!                     clockwise, wrt triangle orientation)
!
!     G1RecTri : septic Bezier triangle
!       Rdata(0:107) - control points (36 x 3)
!
!     CylTri   : image of a linear triangle through a global
!                system of cylindrical coordinates
!
!  VertNo  - vertex points numbers
!
!  EdgeNo  - edge curves numbers with a sign indicating orientation
!     +ncurve : consistency b/w curve (global) orientation and edge
!               local orientation
!     -ncurve : inconstency
!
!  BlockNo - nicknames of blocks adjacent to the figure
!     prism       : nick = nprism*10 + 1
!   ( hexahedron  : nick = nbrick*10 + 2 )
!     tetrahedron : nick = ntetra*10 + 3
!     pyramid     : nick = npyram*10 + 4
!     no block (figure on the boundary) : nick=0
!
!  Domain  - part of a 2D manifold the triangle is on; useful for
!            setting up boundary or interface conditions
!
!  Idata   - integer data associated to the triangle
!
!  Rdata   - real    data associated to the triangle
!
!----------------------------------------------------------------------
      type triangle
        character(len=10)                       :: Type
        integer                                 :: VertNo(3)
        integer                                 :: EdgeNo(3)
        integer                                 :: BlockNo(2)
        integer                                 :: Domain
        integer,          dimension(:), pointer :: Idata
        real(8),          dimension(:), pointer :: Rdata
      endtype triangle
!
!
!----------------------------------------------------------------------
!  RECTANGLE                                                          |
!----------------------------------------------------------------------
!
!  Type =
!
!     BilQua  : bilinear (straight edges) rectangle
!
!     TraQua  : transfinite interpolation rectangle
!
!     PTIRec  : parametric transfinite interpolation rectangle
!       Idata(  1) - conforming surface number
!
!     ImpRec  : implicit rectangle
!       Idata(  1) - conforming surface number
!       Idata(2:5) - bounding surfaces numbers (listed counter
!                     clockwise, wrt rectangle orientation)
!
!     HermRec : geometry reconstruction triangle
!       Rdata - degrees of freedom
!
!     CylRec  : image of a linear rectangle through a global
!               system of cylindrical coordinates
!
!  VertNo  - vertex points numbers
!
!  EdgeNo  - edge curves numbers with a sign indicating orientation
!     +ncurve : consistency b/w curve (global) orientation and edge
!               local orientation
!     -ncurve : inconstency
!
!  BlockNo - nicknames of blocks adjacent to the figure
!     prism       : nick = nprism*10 + 1
!     hexahedron  : nick = nbrick*10 + 2
!   ( tetrahedron : nick = ntetra*10 + 3 )
!     pyramid     : nick = npyram*10 + 4
!     no block (figure on the boundary) : nick=0
!
!  Domain  - part of a 2D manifold the rectangle is on; useful for
!            setting up boundary or interface conditions
!
!  Idata   - integer data associated to the rectangle
!
!  Rdata   - real    data associated to the rectangle
!
!----------------------------------------------------------------------
      type rectangle
        character(len=10)                       :: Type
        integer                                 :: VertNo(4)
        integer                                 :: EdgeNo(4)
        integer                                 :: BlockNo(2)
        integer                                 :: Domain
        integer,          dimension(:), pointer :: Idata
        real(8),          dimension(:), pointer :: Rdata
      endtype rectangle
!
!
!----------------------------------------------------------------------
!  PRISM                                                              |
!----------------------------------------------------------------------
!
!  Type =
!
!     Linear  : trilinear (straight edges) prism
!
!     TIprism : transfinite interpolation prism
!
!  VertNo - vertex points numbers
!
!  EdgeNo - edge curves numbers with a sign indicating orientation
!     +ncurve : consistency b/w curve (global) orientation and edge
!               local orientation
!     -ncurve : inconsistency
!
!  FigNo  - face figures nicknames; the orientations is the figure
!           (global) orientation wrt to the face local orientation
!     horizontal faces : nick = ntrian*10 + norient [0,1,2,3,4,5]
!     vertical   faces : nick = nrecta*10 + norient [0,1,2,3,4,5,6,7]
!
!  Domain - domain number
!
!  Idata  - integer data for the prism
!
!----------------------------------------------------------------------
      type pprism
        character(len=10)              :: Type
        integer                        :: VertNo(6)
        integer                        :: EdgeNo(9)
        integer                        :: FigNo(5)
        integer                        :: Domain
        integer, dimension(:), pointer :: Idata
      endtype pprism
!
!
!----------------------------------------------------------------------
!  HEXAHEDRON                                                         |
!----------------------------------------------------------------------
!
!  Type =
!
!     Linear : trilinear (straight edges) hexahedron
!
!     TraHex : transfinite interpolation hexahedron
!
!     CylHex : image of a linear hex through a global
!              system of cylindrical coordinates
!
!  VertNo - vertex points numbers
!
!  EdgeNo - edge curves numbers with a sign indicating orientation
!     +ncurve : consistency b/w curve (global) orientation and edge
!               local orientation
!     -ncurve : inconsistency
!
!  FigNo  - face figures nicknames; the orientations is the figure
!           (global) orientation wrt to the face local orientation
!     nick = nrecta*10 + norient [0,1,2,3,4,5,6,7]
!
!  Domain - domain number
!
!  Idata  - integer data for the hexahedron
!
!----------------------------------------------------------------------
      type hhexa
        character(len=10)              :: Type
        integer                        :: VertNo(8)
        integer                        :: EdgeNo(12)
        integer                        :: FigNo(6)
        integer                        :: Domain
        integer, dimension(:), pointer :: Idata
      endtype hhexa
!
!
!----------------------------------------------------------------------
!  TETRAHEDRON                                                        |
!----------------------------------------------------------------------
!
!  Type =
!
!     Linear : linear tetrahedron
!
!     TraTet : transfinite interpolation tetrahedron
!
!     CylTet : image of a linear tet through a global
!              system of cylindrical coordinates
!
!  VertNo - vertex points numbers
!
!  EdgeNo - edge curves numbers with a sign indicating orientation
!     +ncurve : consistency b/w curve (global) orientation and edge
!               local orientation
!     -ncurve : inconsistency
!
!  FigNo  - face figures nicknames; the orientations is the figure
!           (global) orientation wrt to the face local orientation
!     nick = ntrian*10 + norient [0,1,2,3,4,5]
!
!  Domain - domain number
!
!  Idata  - integer data for the tetrahedron
!
!----------------------------------------------------------------------
      type ttetra
        character(len=10)              :: Type
        integer                        :: VertNo(4)
        integer                        :: EdgeNo(6)
        integer                        :: FigNo(4)
        integer                        :: Domain
        integer, dimension(:), pointer :: Idata
      endtype ttetra
!
!
!----------------------------------------------------------------------
!  PYRAMID                                                            |
!----------------------------------------------------------------------
!
!  Type =
!
!     Linear  : pyramid with straight edges
!
!     TIPyram : transfinite interpolation pyramid
!
!  VertNo - vertex points numbers
!
!  EdgeNo - edge curves numbers with a sign indicating orientation
!     +ncurve : consistency b/w curve (global) orientation and edge
!               local orientation
!     -ncurve : inconsistency
!
!  FigNo  - face figures nicknames; the orientations is the figure
!           (global) orientation wrt to the face local orientation
!     lateral faces : nick = ntrian*10 + norient [0,1,2,3,4,5]
!     bottom  face  : nick = nrecta*10 + norient [0,1,2,3,4,5,6,7]
!
!  Domain - domain number
!
!  Idata  - integer data for the pyramid
!
!----------------------------------------------------------------------
      type ppyramid
        character(len=10)              :: Type
        integer                        :: VertNo(5)
        integer                        :: EdgeNo(8)
        integer                        :: FigNo(5)
        integer                        :: Domain
        integer, dimension(:), pointer :: Idata
      endtype ppyramid
!
!
!----------------------------------------------------------------------
!  ARRAYS OF 0,1,2,3-D OBJECTS                                        |
!----------------------------------------------------------------------
      type(surface  ), allocatable :: SURFACES(:)
      type(point    ), allocatable :: POINTS(:)
      type(ccurve   ), allocatable :: CURVES(:)
      type(triangle ), allocatable :: TRIANGLES(:)
      type(rectangle), allocatable :: RECTANGLES(:)
      type(hhexa    ), allocatable :: HEXAS(:)
      type(ttetra   ), allocatable :: TETRAS(:)
      type(pprism   ), allocatable :: PRISMS(:)
      type(ppyramid ), allocatable :: PYRAMIDS(:)
!
!
      contains
!
!  ...allocates GMP data structures
   subroutine alloc_GMP
      integer :: i
!
      allocate(SURFACES(MAXSU))
      do i=1,MAXSU
        SURFACES(i)%Type = 'None'
        nullify(SURFACES(i)%Idata)
        nullify(SURFACES(i)%Rdata)
      enddo
!
      allocate(POINTS(MAXNP))
      do i=1,MAXNP
        POINTS(i)%Type = 'None'
        POINTS(i)%NrCurv = 0
        nullify(POINTS(i)%CurvNo)
        nullify(POINTS(i)%Idata)
        nullify(POINTS(i)%Rdata)
      enddo
!
      allocate(CURVES(MAXNC))
      do i=1,MAXNC
        CURVES(i)%Type = 'None'
        CURVES(i)%EndPoNo = 0
        CURVES(i)%NrFig = 0
        nullify(CURVES(i)%FigNo)
        nullify(CURVES(i)%Idata)
        nullify(CURVES(i)%Rdata)
      enddo
!
      allocate(TRIANGLES(MAXTR))
      do i=1,MAXTR
        TRIANGLES(i)%Type = 'None'
        TRIANGLES(i)%VertNo = 0
        TRIANGLES(i)%EdgeNo = 0
        TRIANGLES(i)%BlockNo = 0
        TRIANGLES(i)%Domain = 0
        nullify(TRIANGLES(i)%Idata)
      enddo
!
      allocate(RECTANGLES(MAXRE))
      do i=1,MAXRE
        RECTANGLES(i)%Type = 'None'
        RECTANGLES(i)%VertNo = 0
        RECTANGLES(i)%EdgeNo = 0
        RECTANGLES(i)%BlockNo = 0
        RECTANGLES(i)%Domain = 0
        nullify(RECTANGLES(i)%Idata)
        nullify(RECTANGLES(i)%Rdata)
      enddo
!
      allocate(PRISMS(MAXBT))
      do i=1,MAXBT
        PRISMS(i)%Type = 'None'
        PRISMS(i)%VertNo = 0
        PRISMS(i)%EdgeNo = 0
        PRISMS(i)%FigNo = 0
        PRISMS(i)%Domain = 0
      enddo
!
      allocate(HEXAS(MAXHE))
      do i=1,MAXHE
        HEXAS(i)%Type = 'None'
        HEXAS(i)%VertNo = 0
        HEXAS(i)%EdgeNo = 0
        HEXAS(i)%FigNo = 0
        HEXAS(i)%Domain = 0
        nullify(HEXAS(i)%Idata)
      enddo
!
      allocate(TETRAS(MAXTE))
      do i=1,MAXTE
        TETRAS(i)%Type = 'None'
        TETRAS(i)%VertNo = 0
        TETRAS(i)%EdgeNo = 0
        TETRAS(i)%FigNo = 0
        TETRAS(i)%Domain = 0
        nullify(TETRAS(i)%Idata)
      enddo
!
      allocate(PYRAMIDS(MAXPY))
      do i=1,MAXPY
        PYRAMIDS(i)%Type = 'None'
        PYRAMIDS(i)%VertNo = 0
        PYRAMIDS(i)%EdgeNo = 0
        PYRAMIDS(i)%FigNo = 0
        PYRAMIDS(i)%Domain = 0
        nullify(PYRAMIDS(i)%Idata)
      enddo
!
   end subroutine alloc_GMP
!
!
!----------------------------------------------------------------------
!
!  ...dump out GMP data structure
   subroutine dumpout_GMP
!
      integer, parameter :: ndump=31
!
      integer :: nc,nh,nl,nn,np,nr,ns,nt
      integer :: npri,npyr,ntet
!
      integer :: iprint
      iprint=0
!
      open(unit=ndump,file='files/dumpGMP', &
           form='formatted',access='sequential',status='unknown')
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping GENERAL PARAMETRS'
      write(ndump,*) &
                 NDIM,MANDIM, &
                 NRSURFS, &
                 NRPOINT,NRCURVE,NRTRIAN,NRRECTA, &
                 NRHEXAS,NRTETRA,NRPRISM,NRPYRAM, &
                 NRDOMAIN
!
      write(ndump,*) &
                 MAXSU, &
                 MAXNP,MAXNC,MAXTR,MAXRE, &
                 MAXBT,MAXHE,MAXTE,MAXPY

!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping SURFACES'
      do ns=1,MAXSU
        if (iprint.eq.1) write(*,*) 'dumpGMP: ns = ',ns
        write(ndump,*) SURFACES(ns)%Type
        if (iprint.eq.1) write(*,*) 'dumpGMP: SURFACES(ns)%Type = ', &
                                              SURFACES(ns)%Type
!  .....Idata
        if (associated(SURFACES(ns)%Idata)) then
          nn = ubound(SURFACES(ns)%Idata,1)
          nl = lbound(SURFACES(ns)%Idata,1)
          if (iprint.eq.1) write(*,*) 'dumpGMP: nn FOR Idata= ',nn
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) SURFACES(ns)%Idata
        else
          write(ndump,*) 0
        endif
!  .....Rdata
        if (associated(SURFACES(ns)%Rdata)) then
          nn = ubound(SURFACES(ns)%Rdata,1)
          nl = lbound(SURFACES(ns)%Rdata,1)
          if (iprint.eq.1) write(*,*) 'dumpGMP: nn FOR Rdata= ',nn
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) SURFACES(ns)%Rdata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping POINTS'
      do np=1,MAXNP
        write(ndump,*) POINTS(np)%Type
        write(ndump,*) POINTS(np)%NrCurv
        if (associated(POINTS(np)%CurvNo)) then
          nn = ubound(POINTS(np)%CurvNo,1)
          write(ndump,*) nn
          write(ndump,*) POINTS(np)%CurvNo
        else
          write(ndump,*) 0
        endif
!  .....Idata
        if (associated(POINTS(np)%Idata)) then
          nn = ubound(POINTS(np)%Idata,1)
          nl = lbound(POINTS(np)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) POINTS(np)%Idata
        else
          write(ndump,*) 0
        endif
!  .....Rdata
        if (associated(POINTS(np)%Rdata)) then
          nn = ubound(POINTS(np)%Rdata,1)
          nl = lbound(POINTS(np)%Rdata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) POINTS(np)%Rdata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping CURVES'
      do nc=1,MAXNC
        write(ndump,*) CURVES(nc)%Type
        write(ndump,*) CURVES(nc)%EndPoNo(1:2)
        write(ndump,*) CURVES(nc)%NrFig
        if (associated(CURVES(nc)%FigNo)) then
          nn = ubound(CURVES(nc)%FigNo,1)
          write(ndump,*) nn
          write(ndump,*) CURVES(nc)%FigNo
        else
          write(ndump,*) 0
        endif
!  .....Idata
        if (associated(CURVES(nc)%Idata)) then
          nn = ubound(CURVES(nc)%Idata,1)
          nl = lbound(CURVES(nc)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) CURVES(nc)%Idata
        else
          write(ndump,*) 0
        endif
!  .....Rdata
        if (associated(CURVES(nc)%Rdata)) then
          nn = ubound(CURVES(nc)%Rdata,1)
          nl = lbound(CURVES(nc)%Rdata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) CURVES(nc)%Rdata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping TRIANGLES'
      do nt=1,MAXTR
        write(ndump,*) TRIANGLES(nt)%Type
        write(ndump,*) TRIANGLES(nt)%VertNo(1:3)
        write(ndump,*) TRIANGLES(nt)%EdgeNo(1:3)
        write(ndump,*) TRIANGLES(nt)%BlockNo(1:2)
        write(ndump,*) TRIANGLES(nt)%Domain
!  .....Idata
        if (associated(TRIANGLES(nt)%Idata)) then
          nn = ubound(TRIANGLES(nt)%Idata,1)
          nl = lbound(TRIANGLES(nt)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) TRIANGLES(nt)%Idata
        else
          write(ndump,*) 0
        endif
!  .....Rdata
        if (associated(TRIANGLES(nt)%Rdata)) then
          nn = ubound(TRIANGLES(nt)%Rdata,1)
          nl = lbound(TRIANGLES(nt)%Rdata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) TRIANGLES(nt)%Rdata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping RECTANGLES'
      do nr=1,MAXRE
        write(ndump,*) RECTANGLES(nr)%Type
        write(ndump,*) RECTANGLES(nr)%VertNo(1:4)
        write(ndump,*) RECTANGLES(nr)%EdgeNo(1:4)
        write(ndump,*) RECTANGLES(nr)%BlockNo(1:2)
        write(ndump,*) RECTANGLES(nr)%Domain
!  .....Idata
        if (associated(RECTANGLES(nr)%Idata)) then
          nn = ubound(RECTANGLES(nr)%Idata,1)
          nl = lbound(RECTANGLES(nr)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) RECTANGLES(nr)%Idata
        else
          write(ndump,*) 0
        endif
!  .....Rdata
        if (associated(RECTANGLES(nr)%Rdata)) then
          nn = ubound(RECTANGLES(nr)%Rdata,1)
          nl = lbound(RECTANGLES(nr)%Rdata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) RECTANGLES(nr)%Rdata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping PRISMS'
      do npri=1,MAXBT
        write(ndump,*) PRISMS(npri)%Type
        write(ndump,*) PRISMS(npri)%VertNo(1:6)
        write(ndump,*) PRISMS(npri)%EdgeNo(1:9)
        write(ndump,*) PRISMS(npri)%FigNo(1:5)
        write(ndump,*) PRISMS(npri)%Domain
!  .....Idata
        if (associated(PRISMS(npri)%Idata)) then
          nn = ubound(PRISMS(npri)%Idata,1)
          nl = lbound(PRISMS(npri)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) PRISMS(npri)%Idata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping HEXAS'
      do nh=1,MAXHE
        write(ndump,*) HEXAS(nh)%Type
        write(ndump,*) HEXAS(nh)%VertNo(1:8)
        write(ndump,*) HEXAS(nh)%EdgeNo(1:12)
        write(ndump,*) HEXAS(nh)%FigNo(1:6)
        write(ndump,*) HEXAS(nh)%Domain
!  .....Idata
        if (associated(HEXAS(nh)%Idata)) then
          nn = ubound(HEXAS(nh)%Idata,1)
          nl = lbound(HEXAS(nh)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) HEXAS(nh)%Idata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping TETRAS'
      do ntet=1,MAXTE
        write(ndump,*) TETRAS(ntet)%Type
        write(ndump,*) TETRAS(ntet)%VertNo(1:4)
        write(ndump,*) TETRAS(ntet)%EdgeNo(1:6)
        write(ndump,*) TETRAS(ntet)%FigNo(1:4)
        write(ndump,*) TETRAS(ntet)%Domain
!  .....Idata
        if (associated(TETRAS(ntet)%Idata)) then
          nn = ubound(TETRAS(ntet)%Idata,1)
          nl = lbound(TETRAS(ntet)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) TETRAS(ntet)%Idata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping PYRAMIDS'
      do npyr=1,MAXPY
        write(ndump,*) PYRAMIDS(npyr)%Type
        write(ndump,*) PYRAMIDS(npyr)%VertNo(1:5)
        write(ndump,*) PYRAMIDS(npyr)%EdgeNo(1:8)
        write(ndump,*) PYRAMIDS(npyr)%FigNo(1:5)
        write(ndump,*) PYRAMIDS(npyr)%Domain
!  .....Idata
        if (associated(PYRAMIDS(npyr)%Idata)) then
          nn = ubound(PYRAMIDS(npyr)%Idata,1)
          nl = lbound(PYRAMIDS(npyr)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) PYRAMIDS(npyr)%Idata
        else
          write(ndump,*) 0
        endif
      enddo
!
      if (iprint.eq.1) write(*,*) 'dumpGMP: closing file'
      close(ndump)
!
   end subroutine dumpout_GMP
!
!----------------------------------------------------------------------
!
!  ...dump in GMP data structure
   subroutine dumpin_GMP(Fp)
!
      character(len=*), intent(in)  :: Fp
      integer, parameter :: ndump=31
!
      integer :: nc,nh,nl,nn,np,nr,ns,nt
      integer :: npri,npyr,ntet
!
      integer :: iprint
      iprint=0
!
      open(unit=ndump,file=Fp, &
           form='formatted',access='sequential',status='unknown')
!
      read(ndump,*) &
                 NDIM,MANDIM, &
                 NRSURFS, &
                 NRPOINT,NRCURVE,NRTRIAN,NRRECTA, &
                 NRHEXAS,NRTETRA,NRPRISM,NRPYRAM, &
                 NRDOMAIN
!
      read(ndump,*) &
                 MAXSU, &
                 MAXNP,MAXNC,MAXTR,MAXRE, &
                 MAXBT,MAXHE,MAXTE,MAXPY
!
!!!      if ((MAXSU.ne.MAXSU_old).or. &
!!!          (MAXNP.ne.MAXNP_old).or. &
!!!          (MAXNC.ne.MAXNC_old).or. &
!!!          (MAXTR.ne.MAXTR_old).or. &
!!!          (MAXRE.ne.MAXRE_old).or. &
!!!          (MAXBT.ne.MAXBT_old).or. &
!!!          (MAXHE.ne.MAXHE_old).or. &
!!!          (MAXTE.ne.MAXTE_old).or. &
!!!          (MAXPY.ne.MAXPY_old)) then
!!!        write(*,*) 'dumpin_GMP: INCOMPATIBLE GMP DATA STRUCTURES'
!!!        write(*,*)'MAXSU_old =',MAXSU_old
!!!        write(*,*)'MAXNP_old =',MAXNP_old
!!!        write(*,*)'MAXNC_old =',MAXNC_old
!!!        write(*,*)'MAXTR_old =',MAXTR_old
!!!        write(*,*)'MAXRE_old =',MAXRE_old
!!!        write(*,*)'MAXBT_old =',MAXBT_old
!!!        write(*,*)'MAXHE_old =',MAXHE_old
!!!        write(*,*)'MAXTE_old =',MAXTE_old
!!!        write(*,*)'MAXPY_old =',MAXPY_old
!!!        stop 1
!!!      endif
!
      if (iprint.eq.1) then
        write(*,*)'---------------------------------------------'
        write(*,*)'dumpin_GMP: '
        write(*,8000)NDIM
8000    format('     NDIM     = ',i12)
        write(*,8001)MANDIM
8001    format('     MANDIM   = ',i12)
        write(*,8002)NRSURFS
8002    format('     NRSURFS  = ',i12)
        write(*,8003)NRPOINT
8003    format('     NRPOINT  = ',i12)
        write(*,8004)NRCURVE
8004    format('     NRCURVE  = ',i12)
        write(*,8005)NRTRIAN
8005    format('     NRTRIAN  = ',i12)
        write(*,8006)NRRECTA
8006    format('     NRRECTA  = ',i12)
        write(*,8007)NRHEXAS
8007    format('     NRHEXAS  = ',i12)
        write(*,8008)NRTETRA
8008    format('     NRTETRA  = ',i12)
        write(*,8009)NRPRISM
8009    format('     NRPRIAM  = ',i12)
        write(*,8010)NRPYRAM
8010    format('     NRPYRAM  = ',i12)
        write(*,8011)NRDOMAIN
8011    format('     NRDOMAIN = ',i12)
        write(*,*)'---------------------------------------------'
        call pause
      endif
!
      call alloc_GMP
!
      do ns=1,MAXSU
        read(ndump,*) SURFACES(ns)%Type
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(SURFACES(ns)%Idata(nl:nn))
          read(ndump,*) SURFACES(ns)%Idata
        else
          nullify(SURFACES(ns)%Idata)
        endif
!  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(SURFACES(ns)%Rdata(nl:nn))
          read(ndump,*) SURFACES(ns)%Rdata
        else
          nullify(SURFACES(ns)%Rdata)
        endif
      enddo
!
      do np=1,MAXNP
        read(ndump,*) POINTS(np)%Type
        read(ndump,*) POINTS(np)%NrCurv
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(POINTS(np)%CurvNo(nn))
          read(ndump,*) POINTS(np)%CurvNo
        else
          nullify(POINTS(np)%CurvNo)
        endif
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(POINTS(np)%Idata(nl:nn))
          read(ndump,*) POINTS(np)%Idata
        else
          nullify(POINTS(np)%Idata)
        endif
!  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(POINTS(np)%Rdata(nl:nn))
          read(ndump,*) POINTS(np)%Rdata
        else
          nullify(POINTS(np)%Rdata)
        endif
      enddo
!
      do nc=1,MAXNC
        read(ndump,*) CURVES(nc)%Type
        read(ndump,*) CURVES(nc)%EndPoNo(1:2)
        read(ndump,*) CURVES(nc)%NrFig
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(CURVES(nc)%FigNo(nn))
          read(ndump,*) CURVES(nc)%FigNo
        else
          nullify(CURVES(nc)%FigNo)
        endif
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(CURVES(nc)%Idata(nl:nn))
          read(ndump,*) CURVES(nc)%Idata
        else
          nullify(CURVES(nc)%Idata)
        endif
!  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(CURVES(nc)%Rdata(nl:nn))
          read(ndump,*) CURVES(nc)%Rdata
        else
          nullify(CURVES(nc)%Rdata)
        endif
      enddo
!
      do nt=1,MAXTR
        read(ndump,*) TRIANGLES(nt)%Type
        read(ndump,*) TRIANGLES(nt)%VertNo(1:3)
        read(ndump,*) TRIANGLES(nt)%EdgeNo(1:3)
        read(ndump,*) TRIANGLES(nt)%BlockNo(1:2)
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(TRIANGLES(nt)%Idata(nl:nn))
          read(ndump,*) TRIANGLES(nt)%Idata
        else
          nullify(TRIANGLES(nt)%Idata)
        endif
!  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(TRIANGLES(nt)%Rdata(nl:nn))
          read(ndump,*) TRIANGLES(nt)%Rdata
        else
          nullify(TRIANGLES(nt)%Rdata)
        endif
      enddo
!
      do nr=1,MAXRE
        read(ndump,*) RECTANGLES(nr)%Type
        read(ndump,*) RECTANGLES(nr)%VertNo(1:4)
        read(ndump,*) RECTANGLES(nr)%EdgeNo(1:4)
        read(ndump,*) RECTANGLES(nr)%BlockNo(1:2)
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(RECTANGLES(nr)%Idata(nl:nn))
          read(ndump,*) RECTANGLES(nr)%Idata
        else
          nullify(RECTANGLES(nr)%Idata)
        endif
!  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(RECTANGLES(nr)%Rdata(nl:nn))
          read(ndump,*) RECTANGLES(nr)%Rdata
        else
          nullify(RECTANGLES(nr)%Rdata)
        endif
      enddo
!
      do npri=1,MAXBT
        read(ndump,*) PRISMS(npri)%Type
        read(ndump,*) PRISMS(npri)%VertNo(1:6)
        read(ndump,*) PRISMS(npri)%EdgeNo(1:9)
        read(ndump,*) PRISMS(npri)%FigNo(1:5)
        read(ndump,*) PRISMS(npri)%Domain
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nn
          allocate(PRISMS(npri)%Idata(nl:nn))
          read(ndump,*) PRISMS(npri)%Idata
        else
          nullify(PRISMS(npri)%Idata)
        endif
      enddo
!
      do nh=1,MAXHE
        read(ndump,*) HEXAS(nh)%Type
        read(ndump,*) HEXAS(nh)%VertNo(1:8)
        read(ndump,*) HEXAS(nh)%EdgeNo(1:12)
        read(ndump,*) HEXAS(nh)%FigNo(1:6)
        read(ndump,*) HEXAS(nh)%Domain
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(HEXAS(nh)%Idata(nl:nn))
          read(ndump,*) HEXAS(nh)%Idata
        else
          nullify(HEXAS(nh)%Idata)
        endif
      enddo
!
      do ntet=1,MAXTE
        read(ndump,*) TETRAS(ntet)%Type
        read(ndump,*) TETRAS(ntet)%VertNo(1:4)
        read(ndump,*) TETRAS(ntet)%EdgeNo(1:6)
        read(ndump,*) TETRAS(ntet)%FigNo(1:4)
        read(ndump,*) TETRAS(ntet)%Domain
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(TETRAS(ntet)%Idata(nl:nn))
          read(ndump,*) TETRAS(ntet)%Idata
        else
          nullify(TETRAS(ntet)%Idata)
        endif
      enddo
!
      do npyr=1,MAXPY
        read(ndump,*) PYRAMIDS(npyr)%Type
        read(ndump,*) PYRAMIDS(npyr)%VertNo(1:5)
        read(ndump,*) PYRAMIDS(npyr)%EdgeNo(1:8)
        read(ndump,*) PYRAMIDS(npyr)%FigNo(1:5)
        read(ndump,*) PYRAMIDS(npyr)%Domain
!  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(PYRAMIDS(npyr)%Idata(nl:nn))
          read(ndump,*) PYRAMIDS(npyr)%Idata
        else
          nullify(PYRAMIDS(npyr)%Idata)
        endif
      enddo
!
      close(ndump)
!
   end subroutine dumpin_GMP
!
!----------------------------------------------------------------------
!
   subroutine print_GMP
!
      use control
      implicit none
!
      integer :: i,idec,nc,nh,np,nr,ns,nt
      integer :: npri,npyr,ntet,nrcurv,nrfig
!
 10   continue
      write(*,*) 'print_GMP: SELECT:'
      write(*,*) 'EXIT....................................0'
      write(*,*) 'GENERAL DATA STRUCTURE PARAMETERS.......1'
      write(*,*) 'POINT DATA..............................2'
      write(*,*) 'CURVE DATA..............................3'
      write(*,*) 'TRIANGLE DATA...........................4'
      write(*,*) 'RECTANGLE DATA..........................5'
      write(*,*) 'PRISM DATA..............................6'
      write(*,*) 'HEXAHEDRON DATA.........................7'
      write(*,*) 'TETRAHEDRON DATA........................8'
      write(*,*) 'PYRAMID DATA............................9'
      write(*,*) 'SURFACE DATA...........................10'
      read(*,*) idec
!
      select case(idec)
!
      case(0)
        return
!
      case(1)
        write(*,7001) NDIM,MANDIM
 7001   format(' NDIM,MANDIM = ',2i3)
        write(*,7002) NRSURFS
 7002   format(' NRSURFS = ',i10)
        write(*,7003) NRPOINT,NRCURVE,NRTRIAN,NRRECTA
 7003   format(' NRPOINT,NRCURVE,NRTRIAN,NRRECTA = ',4i8)
        write(*,7004) NRHEXAS,NRTETRA,NRPRISM,NRPYRAM
 7004   format(' NRHEXAS,NRTETRA,NRPRISM,NRPYRAM = ',4i8)
        write(*,7005) NRDOMAIN
 7005   format(' NRDOMAIN = ',i3)
        write(*,7006) MAXSU, &
                      MAXNP,MAXNC,MAXTR,MAXRE, &
                      MAXBT,MAXHE,MAXTE,MAXPY
 7006   format(' MAXSU,MAXNP,MAXNC,MAXTR,MAXRE,MAXBT,MAXHE,MAXTE,', &
               'MAXPY = ',9i7)
!
      case(2)
        write(*,*) 'SET POINT NUMBER'
        read(*,*) np
        if ((np.le.0).or.(np.gt.NRPOINT)) goto 10
        write(*,7101) np, POINTS(np)%Type
 7101   format(' POINT = ',i8,' TYPE = ',a10)
        nrcurv = POINTS(np)%NrCurv
        write(*,7102) POINTS(np)%CurvNo(1:nrcurv)
 7102   format(' CONNECTED CURVES = ',10i8)
        if (associated(POINTS(np)%Idata)) then
          write(*,8002) POINTS(np)%Idata
        endif
        if (associated(POINTS(np)%Rdata)) then
          write(*,8003) POINTS(np)%Rdata
        endif
!
      case(3)
        write(*,*) 'SET CURVE NUMBER'
        read(*,*) nc
        if ((nc.le.0).or.(nc.gt.NRCURVE)) goto 10
        write(*,7201) nc, CURVES(nc)%Type
 7201   format(' CURVE = ',i8,' TYPE = ',a10)
        write(*,7202) CURVES(nc)%EndPoNo(1:2)
 7202   format(' END POINTS = ',2i8)
        nrfig = CURVES(nc)%NrFig
        write(*,7203) CURVES(nc)%FigNo(1:nrfig)
 7203   format(' CONNECTED FIGURES = ',10i8)
        if (associated(CURVES(nc)%Idata)) then
          write(*,8002) CURVES(nc)%Idata
        endif
        if (associated(CURVES(nc)%Rdata)) then
          select case(CURVES(nc)%Type)
          case('QuaCir','QuaEl1','QuaEl2')
            write(*,5555) CURVES(nc)%Rdata(1:3)
 5555       format(' Rdata(1:3) (center) = ',3(e12.5,2x))
          case('QuaSEl')
            write(*,5556) CURVES(nc)%Rdata(1:3)
 5556       format(' Rdata(1:3) (center) = ',3(e12.5,2x))
            write(*,5557) CURVES(nc)%Rdata(4:5)
 5557       format(' Rdata(4:5) (nx,ny) = ',2(e12.5,2x))
          case default
            do i=0,7
              write(*,7405)i,CURVES(nc)%Rdata(3*i:3*i+2)
            enddo
          endselect
        endif
!
      case(4)
        write(*,*) 'SET TRIANGLE NUMBER'
        read(*,*) nt
        if ((nt.le.0).or.(nt.gt.NRTRIAN)) goto 10
        write(*,7501) nt,TRIANGLES(nt)%Type
 7501   format(' TRIANGLE = ',i8,' TYPE = ',a10)
        write(*,7402) TRIANGLES(nt)%VertNo(1:3)
        write(*,7403) TRIANGLES(nt)%EdgeNo(1:3)
        write(*,7404) TRIANGLES(nt)%BlockNo(1:2)
        if (associated(TRIANGLES(nt)%Idata)) then
          write(*,8002) TRIANGLES(nt)%Idata
        endif
        if (associated(TRIANGLES(nt)%Rdata)) then
          if (TRIANGLES(nt)%Type.eq.'G1RecTri') then
          write(*,*)'CONTROL POINTS'
            do i = 0,35
              write(*,7405)i,TRIANGLES(nt)%Rdata(3*i:3*i+2)
 7405         format(' b',i2,' = ',3(e20.13,2x))
            enddo
          endif
        endif
!
      case(5)
        write(*,*) 'SET RECTANGLE NUMBER'
        read(*,*) nr
        if ((nr.le.0).or.(nr.gt.NRRECTA)) goto 10
        write(*,7401) nr,RECTANGLES(nr)%Type
 7401   format(' RECTANGLE = ',i8,' TYPE = ',a10)
        write(*,7402) RECTANGLES(nr)%VertNo(1:4)
 7402   format('VERTICES = ',4i8)
        write(*,7403) RECTANGLES(nr)%EdgeNo(1:4)
 7403   format('EDGES = ',4i8)
        write(*,7404) RECTANGLES(nr)%BlockNo(1:2)
 7404   format(' ADJACENT BLOCKS = ',2i8)
        if (associated(RECTANGLES(nr)%Idata)) then
          write(*,8002) RECTANGLES(nr)%Idata
        endif
        if (associated(RECTANGLES(nr)%Rdata)) then
          write(*,8003) RECTANGLES(nr)%Rdata
        endif
!
      case(6)
        write(*,*) 'SET PRISM NUMBER'
        read(*,*) npri
        if ((npri.le.0).or.(npri.gt.NRPRISM)) goto 10
        write(*,7701) npri,PRISMS(npri)%Type
 7701   format(' PRISM = ',i8,' TYPE = ',a10)
        write(*,7602) PRISMS(npri)%VertNo(1:6)
        write(*,7603) PRISMS(npri)%EdgeNo(1:9)
        write(*,7604) PRISMS(npri)%FigNo(1:5)
        write(*,7605) PRISMS(npri)%Domain
!
      case(7)
        write(*,*) 'SET HEXAHEDRON NUMBER'
        read(*,*) nh
        if ((nh.le.0).or.(nh.gt.NRHEXAS)) goto 10
        write(*,7601) nh,HEXAS(nh)%Type
 7601   format(' HEXA = ',i8,' TYPE = ',a10)
        write(*,7602) HEXAS(nh)%VertNo(1:8)
 7602   format(' VERTICES = ',8i8)
        write(*,7603) HEXAS(nh)%EdgeNo(1:12)
 7603   format(' EDGES = ',12i8)
        write(*,7604) HEXAS(nh)%FigNo(1:6)
 7604   format(' FACES = ',6i8)
        write(*,7605) HEXAS(nh)%Domain
 7605   format(' Domain = ',i8)
       if (associated(HEXAS(nh)%Idata)) then
          write(*,8002) HEXAS(nh)%Idata
        endif
!
      case(8)
        write(*,*) 'SET TETRAHEDRON NUMBER'
        read(*,*) ntet
        if ((ntet.le.0).or.(ntet.gt.NRTETRA)) goto 10
        write(*,7801) ntet,TETRAS(ntet)%Type
 7801   format(' TETRA = ',i8,' TYPE = ',a10)
        write(*,7602) TETRAS(ntet)%VertNo(1:4)
        write(*,7603) TETRAS(ntet)%EdgeNo(1:6)
        write(*,7604) TETRAS(ntet)%FigNo(1:4)
        write(*,7605) TETRAS(ntet)%Domain
        if (associated(TETRAS(ntet)%Idata)) then
          write(*,8002) TETRAS(ntet)%Idata
        endif
!
      case(9)
        write(*,*) 'SET PYRAMID NUMBER'
        read(*,*) npyr
        if ((npyr.le.0).or.(npyr.gt.NRPYRAM)) goto 10
        write(*,7901) npyr,PYRAMIDS(npyr)%Type
 7901   format(' PYRAMID = ',i8,' TYPE = ',a10)
        write(*,7602) PYRAMIDS(npyr)%VertNo(1:5)
        write(*,7603) PYRAMIDS(npyr)%EdgeNo(1:8)
        write(*,7604) PYRAMIDS(npyr)%FigNo(1:5)
        write(*,7605) PYRAMIDS(npyr)%Domain
        if (associated(PYRAMIDS(npyr)%Idata)) then
          write(*,8002) PYRAMIDS(npyr)%Idata
        endif
!
      case(10)
        write(*,*) 'SET SURFACES NUMBER'
        read(*,*) ns
        if ((ns.le.0).or.(ns.gt.NRSURFS)) goto 10
        write(*,8001) ns,SURFACES(ns)%Type
 8001   format(' SURFACE = ',i8,' TYPE = ',a10)
        if (associated(SURFACES(ns)%Idata)) then
          write(*,8002) SURFACES(ns)%Idata
 8002     format(' Idata = ',10i10)
        endif
        if (associated(SURFACES(ns)%Rdata)) then
          write(*,8003) SURFACES(ns)%Rdata
 8003     format(' Rdata = ',8e12.5)
        endif
!
      case default
        goto 10
      end select
      goto 10
!
   end subroutine print_GMP
!
!----------------------------------------------------------------------
!
   subroutine set_gmp_parameters(NDIM_loc,   &
                                    MANDIM_loc, &
                                    MAXSU_loc,  &
                                    MAXNP_loc,  &
                                    MAXNC_loc,  &
                                    MAXTR_loc,  &
                                    MAXRE_loc,  &
                                    MAXBT_loc,  &
                                    MAXHE_loc,  &
                                    MAXTE_loc,  &
                                    MAXPY_loc)
      implicit none
      integer, intent(in) :: NDIM_loc
      integer, intent(in) :: MANDIM_loc
      integer, intent(in) :: MAXSU_loc
      integer, intent(in) :: MAXNP_loc
      integer, intent(in) :: MAXNC_loc
      integer, intent(in) :: MAXTR_loc
      integer, intent(in) :: MAXRE_loc
      integer, intent(in) :: MAXBT_loc
      integer, intent(in) :: MAXHE_loc
      integer, intent(in) :: MAXTE_loc
      integer, intent(in) :: MAXPY_loc
!
      integer, parameter :: iprint=0
!
!  ...GMP control parameters...........................................
      NDIM   = NDIM_loc
      MANDIM = MANDIM_loc
      MAXSU  = MAXSU_loc
      MAXNP  = MAXNP_loc
      MAXNC  = MAXNC_loc
      MAXTR  = MAXTR_loc
      MAXRE  = MAXRE_loc
      MAXBT  = MAXBT_loc
      MAXHE  = MAXHE_loc
      MAXTE  = MAXTE_loc
      MAXPY  = MAXPY_loc
!
!  ...fine grid parameters.............................................
      NRPO_FG  = 0
      MAXNP_FG = MAXNP_loc
!
      if (iprint.eq.1) then
        write(*,*)'-- GMP Control Parameters --'
        write(*,8050)NDIM
 8050   format(' NDIM     = ',i12)
        write(*,8051)MANDIM
 8051   format(' MANDIM   = ',i12)
        write(*,8052)MAXSU
 8052   format(' MAXSU    = ',i12)
        write(*,8053)MAXNP
 8053   format(' NMAXNP   = ',i12)
        write(*,8054)MAXNC
 8054   format(' MAXNC    = ',i12)
        write(*,8055)MAXTR
 8055   format(' MAXTR    = ',i12)
        write(*,8056)MAXRE
 8056   format(' MAXRE    = ',i12)
        write(*,8057)MAXBT
 8057   format(' MAXBT    = ',i12)
        write(*,8058)MAXHE
 8058   format(' MAXHE    = ',i12)
        write(*,8059)MAXTE
 8059   format(' MAXTE    = ',i12)
        write(*,8060)MAXPY
 8060   format(' MAXPY    = ',i12)
        write(*,8061)NRPO_FG
 8061   format(' NRPO_FG  = ',i12)
        write(*,8062)MAXNP_FG
 8062   format(' MAXNP_FG = ',i12)
        write(*,*)''
      endif
!
   end subroutine set_gmp_parameters
!
!----------------------------------------------------------------------
!
   subroutine print_GMP_parameters
      use environment , only : QUIET_MODE
!
      if (.not. QUIET_MODE) then
        write(*,*)'-- GMP Parameters --'
        write(*,8050)NDIM,MANDIM
 8050   format(' NDIM,MANDIM   = ',2(i9,2x))
        write(*,8052)NRSURFS,MAXSU
 8052   format(' NRSURFS,MAXSU = ',2(i9,2x))
        write(*,8053)NRPOINT,MAXNP
 8053   format(' NRPOINT,MAXNP = ',2(i9,2x))
        write(*,8054)NRCURVE,MAXNC
 8054   format(' NRCURVE,MAXNC = ',2(i9,2x))
        write(*,8055)NRTRIAN,MAXTR
 8055   format(' NRTRIAN,MAXTR = ',2(i9,2x))
        write(*,8056)NRRECTA,MAXRE
 8056   format(' NRRECTA,MAXRE = ',2(i9,2x))
        write(*,8057)NRPRISM,MAXBT
 8057   format(' NRPRISM,MAXBT = ',2(i9,2x))
        write(*,8058)NRHEXAS,MAXHE
 8058   format(' NRHEXAS,MAXHE = ',2(i9,2x))
        write(*,8059)NRTETRA,MAXTE
 8059   format(' NRTETRA,MAXTE = ',2(i9,2x))
        write(*,8060)NRPYRAM,MAXPY
 8060   format(' NRPYRAM,MAXPY = ',2(i9,2x))
!!!        write(*,8061)NRPO_FG,MAXNP_FG
!!! 8061   format(' NRPO_FG  = ',i12)
        write(*,*)''
      endif
!
   end subroutine print_GMP_parameters
!
!----------------------------------------------------------------------
!
      integer function Ndomain_GMP(IGMP_block)
!
      integer IGMP_block, no, lab
      call decode(IGMP_block, no,lab)
      select case(lab)
      case(1)
        Ndomain_GMP = PRISMS(no)%Domain
      case(2)
        Ndomain_GMP = HEXAS(no)%Domain
      case(3)
        Ndomain_GMP = TETRAS(no)%Domain
      case(4)
        Ndomain_GMP = PYRAMIDS(no)%Domain
      case default
        write(*,*) 'Ndomain_GMP: lab = ',lab
        stop
      endselect
!
      end function Ndomain_GMP
!
   end module GMP
