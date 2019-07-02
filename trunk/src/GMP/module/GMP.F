c----------------------------------------------------------------------
c> Purporse : data structure arrays for Geometry Modeling Package
c
c  @revision Nov 12
c----------------------------------------------------------------------
c
      module GMP
c
c----------------------------------------------------------------------
c  PARAMETERS                                                         |
c----------------------------------------------------------------------
c     dimension of the problem (2, 3)
      integer, save :: NDIM
c     dimension of the manifold (2, 3)
      integer, save :: MANDIM
c     number of algebraic surfaces
      integer, save :: NRSURFS
c     number of points
      integer, save :: NRPOINT
c     number of curves
      integer, save :: NRCURVE
c     number of triangles
      integer, save :: NRTRIAN
c     number of rectangles
      integer, save :: NRRECTA
c     number of hexahedra
      integer, save :: NRHEXAS
c     number of tetrahedra
      integer, save :: NRTETRA
c     number of prisms
      integer, save :: NRPRISM
c     number of pyramids
      integer, save :: NRPYRAM
c     number of subdomains
      integer, save :: NRDOMAIN
c     number of fine grid points
      integer, save :: NRPO_FG
c     fine grid scale
      double precision, save :: SCALE_FG
c
c     maximum number of algebraic surfaces
      integer, save :: MAXSU
c     maximum number of points
      integer, save :: MAXNP
c     maximum number of curves
      integer, save :: MAXNC
c     maximum number of triangles
      integer, save :: MAXTR
c     maximum number of rectangles
      integer, save :: MAXRE
c     maximum number of prismatic blocks
      integer, save :: MAXBT
c     maximum number of hexahedral blocks
      integer, save :: MAXHE
c     maximum number of tetrahedral blocks
      integer, save :: MAXTE
c     maximum number of pyramids
      integer, save :: MAXPY
c     maximum number of fine grid points
      integer, save :: MAXNP_FG
c
c
c----------------------------------------------------------------------
c  SURFACE                                                            |
c----------------------------------------------------------------------
c
c  Type =
c
c     VecPt     : plane defined by a point and a normal
c       Rdata(1:3) - coordinates of the point
c       Rdata(4:6) - normal
c
c     ThrPt     : plane defined by the coordinates of three points
c       Rdata(1:3) - coordinates of the 1st point
c       Rdata(4:6) - coordinates of the 2nd point
c       Rdata(7:9) - coordinates of the 3rd point
c
c     PtNoVec   : plane defined by a point number and a normal
c       Idata(1)   - point number
c       Rdata(1:3) - normal
c
c     PtNo2Pt   : plane defined by 1 point number and the coordinates
c                 of 2 points
c       Idata(1)   - point number
c       Rdata(1:3) - coordinates of the 1st point
c       Rdata(4:6) - coordinates of the 2nd point
c
c     PPwCC : yz plane parametrized with cylindrical coordinates
c       Rdata(1)   - x coordinate
c       Rdata(2)   - rmin
c       Rdata(3)   - rmax
c
c     Sphere    : sphere
c       Rdata(1:3) - coordinates of the center
c       Rdata(4)   - radius
c
c     Cylinder  : infinite cylinder
c       Rdata(1:3) - coordinates of a point on the axis
c       Rdata(4:6) - axis
c       Rdata(7)   - radius
c
c     Ellipsoid : ellipsoid
c       Rdata(1:3) - coordinates of the center
c       Rdata(4:6) - length of 1st,2nd,3rd semiaxis
c
c     Cone      : cone
c       Rdata(1:3) - coordinates of the vertex
c       Rdata(4:6) - axis
c       Rdata(7)   - half-aperture [Radiants]
c
c     RecSurf   : reconstructed surface
c
c  Idata  - integer data associated to the surface
c
c  Rdata  - real    data associated to the surface
c
c----------------------------------------------------------------------
      type surface
        character(len=10)                       :: Type
        integer,          dimension(:), pointer :: Idata
        double precision, dimension(:), pointer :: Rdata
      endtype surface
c
c
c----------------------------------------------------------------------
c  POINT                                                              |
c----------------------------------------------------------------------
c
c  Type =
c
c     Regular  : a point
c       Rdata(1:3) - coordinates
c
c     Implicit : implicitly defined point
c       Idata(1:3) - surface nubers identifing the point
c       Rdata(1:3) - coordinates determined by NR iterations
c
c     CoorNrm  : point with a normal
c       Rdata(1:3) - coordinates
c       Rdata(1:3) - normal
c
c     SharpPt  : point on a sharp edge
c       Rdata(1:3) - coordinates
c       Rdata(3:6) - normal to the 1st surface defining the sharp edge
c       Rdata(7:9) - normal to the 2nd surface defining the sharp edge
c
c  NrCurv - number of curves connected to the point
c
c  CurvNo - list of connected curves
c
c  Idata  - integer data associated to the point
c
c  Rdata  - real    data associated to the point
c
c----------------------------------------------------------------------
      type point
        character(len=10)                       :: Type
        integer                                 :: NrCurv
        integer,          dimension(:), pointer :: CurvNo
        integer,          dimension(:), pointer :: Idata
        double precision, dimension(:), pointer :: Rdata
      endtype point
c
c
c----------------------------------------------------------------------
c  CURVE                                                              |
c----------------------------------------------------------------------
c
c  Type =
c
c     Seglin    : straight segment
c
c     QuaCir    : quarter of a circle
c       Rdata(1:3) - coordinates of the center
c
c     QuaEl1    : quarter of a ellipse (0 to pi/2)
c       Rdata(1:3) - coordinates of the center
c
c     QuaEl2    : quarter of a ellipse (-pi/4 to pi/4)
c       Rdata(1:3) - coordinates of the center
c
c     QuaSEl    : quarter of a superellipse (0 to pi/2)
c       Rdata(1:3) - coordinates of the center
c       Rdata(4:5) - px,py powers of superellipse
c
c     ImpCir    : implicit curve
c       Idata(1:4) - surfaces defining the curve
c
c     3HermCur  : cubic Hermite curve
c       Rdata(1:3) - derivative at 1st end point
c       Rdata(4:6) - derivative at 2nd end point
c
c     1SurfsCur : curve lying on 1 algebraic surface
c       Idata(1) - surface number
c
c     2SurfsCur : curve on the intersection of 2 algebraic surfaces
c       Idata(1:2) - surface numbers
c
c     3SurfsCur : curve on the intersection of 3 algebraic surfaces
c       Idata(1:3) - surface numbers
c
c     HermCur   : quintic Hermite curve
c       Rdata( 1: 3) - derivative at 1st end point
c       Rdata( 4: 6) - derivative at 2nd end point
c       Rdata( 7: 9) - second derivative at 1st end point
c       Rdata(10:11) - second derivative at 2nd end point
c
c     5Bezier   : quintic Bezier curve
c       Rdata(0:17) - control points (6 x 3) from 1st to 2nd endpoint
c
c     7Bezier   : septic Bezier curve
c       Rdata(0:23) - control points (8 x 3) from 1st to 2nd endpoint
c
c     CylCur    : image of a straight line segment through a global
c                 system of cylindrical coordinates
c
c  EndPoNo - curve endpoints numbers
c
c  NrFig   - number of figures connected to the curve
c
c  FigNo   - list of nicknames of connected figures
c     triangle  : nick = ntrian*10 + 1
c     rectangle : nick = nrecta*10 + 2
c
c  Idata  - integer data associated to the curve
c
c  Rdata  - real    data associated to the curve
c
c----------------------------------------------------------------------
      type ccurve
        character(len=10)                       :: Type
        integer                                 :: EndPoNo(2)
        integer                                 :: NrFig
        integer,          dimension(:), pointer :: FigNo
        integer,          dimension(:), pointer :: Idata
        double precision, dimension(:), pointer :: Rdata
      endtype ccurve
c
c
c----------------------------------------------------------------------
c  TRIANGLE                                                           |
c----------------------------------------------------------------------
c
c  Type =
c
c     PlaneTri : plane triangle
c
c     TransTri : transfinite interpolation triangle
c
c     PTITri   : parametric transfinite interpolation triangle
c       Idata(  1) - conforming surface number
c
c     ImpliTri : implicit triangle
c       Idata(  1) - conforming surface number
c       Idata(2:4) - bounding surfaces numbers (listed counter
c                     clockwise, wrt triangle orientation)
c
c     G1RecTri : septic Bezier triangle
c       Rdata(0:107) - control points (36 x 3)
c
c     SpherTri : LEGACY
c       Idata,Rdata - UNKNOWN
c
c     QtCirTri : LEGACY
c       Idata,Rdata - UNKNOWN
c
c     PaSphTri : LEGACY
c       Idata,Rdata - UNKNOWN
c
c     CylTri    : image of a linear triangle through a global
c                 system of cylindrical coordinates
c
c  VertNo  - vertex points numbers
c
c  EdgeNo  - edge curves numbers with a sign indicating orientation
c     +ncurve : consistency b/w curve (global) orientation and edge
c               local orientation
c     -ncurve : inconstency
c
c  BlockNo - nicknames of blocks adjacent to the figure
c     prism       : nick = nprism*10 + 1
c   ( hexahedron  : nick = nbrick*10 + 2 )
c     tetrahedron : nick = ntetra*10 + 3
c     pyramid     : nick = npyram*10 + 4
c     no block (figure on the boundary) : nick=0
c
c  Idata   - integer data associated to the triangle
c
c  Rdata   - real    data associated to the triangle
c
c----------------------------------------------------------------------
      type triangle
        character(len=10)                       :: Type
        integer                                 :: VertNo(3)
        integer                                 :: EdgeNo(3)
        integer                                 :: BlockNo(2)
        integer,          dimension(:), pointer :: Idata
        double precision, dimension(:), pointer :: Rdata
      endtype triangle
c
c
c----------------------------------------------------------------------
c  RECTANGLE                                                          |
c----------------------------------------------------------------------
c
c  Type =
c
c     BilQua  : bilinear (straight edges) rectangle
c
c     TraQua  : transfinite interpolation rectangle
c
c     PTIRec  : parametric transfinite interpolation rectangle
c       Idata(  1) - conforming surface number
c
c     ImpRec  : implicit rectangle
c       Idata(  1) - conforming surface number
c       Idata(2:5) - bounding surfaces numbers (listed counter
c                     clockwise, wrt rectangle orientation)
c
c     HermRec : geometry reconstruction triangle
c       Rdata - degrees of freedom
c
c
c     CylRec    : image of a linear rectangle through a global
c                 system of cylindrical coordinates
c
c     CylRec (other) : LEGACY (this also exists as a routine written by Paolo)
c       Idata,Rdata - UNKNOWN
c
c     ConfRec : LEGACY
c       Idata,Rdata - UNKNOWN
c
c  VertNo  - vertex points numbers
c
c  EdgeNo  - edge curves numbers with a sign indicating orientation
c     +ncurve : consistency b/w curve (global) orientation and edge
c               local orientation
c     -ncurve : inconstency
c
c  BlockNo - nicknames of blocks adjacent to the figure
c     prism       : nick = nprism*10 + 1
c     hexahedron  : nick = nbrick*10 + 2
c   ( tetrahedron : nick = ntetra*10 + 3 )
c     pyramid     : nick = npyram*10 + 4
c     no block (figure on the boundary) : nick=0
c
c  Idata   - integer data associated to the rectangle
c
c  Rdata   - real    data associated to the rectangle
c
c----------------------------------------------------------------------
      type rectangle
        character(len=10)                       :: Type
        integer                                 :: VertNo(4)
        integer                                 :: EdgeNo(4)
        integer                                 :: BlockNo(2)
        integer,          dimension(:), pointer :: Idata
        double precision, dimension(:), pointer :: Rdata
      endtype rectangle
c
c
c----------------------------------------------------------------------
c  PRISM                                                              |
c----------------------------------------------------------------------
c
c  Type =
c
c     Linear  : trilinear (straight edges) prism
c
c     TIprism : transfinite interpolation prism
c
c  VertNo - vertex points numbers
c
c  EdgeNo - edge curves numbers with a sign indicating orientation
c     +ncurve : consistency b/w curve (global) orientation and edge
c               local orientation
c     -ncurve : inconsistency
c
c  FigNo  - face figures nicknames; the orientations is the figure
c           (global) orientation wrt to the face local orientation
c     horizontal faces : nick = ntrian*10 + norient [0,1,2,3,4,5]
c     vertical   faces : nick = nrecta*10 + norient [0,1,2,3,4,5,6,7]
c
c  Domain - domain number
c
c  Idata  - integer data for the prism
c
c----------------------------------------------------------------------
      type pprism
        character(len=10)              :: Type
        integer                        :: VertNo(6)
        integer                        :: EdgeNo(9)
        integer                        :: FigNo(5)
        integer                        :: Domain
        integer, dimension(:), pointer :: Idata
      endtype pprism
c
c
c----------------------------------------------------------------------
c  HEXAHEDRON                                                         |
c----------------------------------------------------------------------
c
c  Type =
c
c     Linear : trilinear (straight edges) hexahedron
c
c     TraHex : transfinite interpolation hexahedron
c
c     CylHex : image of a linear hex through a global
c              system of cylindrical coordinates
c
c  VertNo - vertex points numbers
c
c  EdgeNo - edge curves numbers with a sign indicating orientation
c     +ncurve : consistency b/w curve (global) orientation and edge
c               local orientation
c     -ncurve : inconsistency
c
c  FigNo  - face figures nicknames; the orientations is the figure
c           (global) orientation wrt to the face local orientation
c     nick = nrecta*10 + norient [0,1,2,3,4,5,6,7]
c
c  Domain - domain number
c
c  Idata  - integer data for the hexahedron
c
c----------------------------------------------------------------------
      type hhexa
        character(len=10)              :: Type
        integer                        :: VertNo(8)
        integer                        :: EdgeNo(12)
        integer                        :: FigNo(6)
        integer                        :: Domain
        integer, dimension(:), pointer :: Idata
      endtype hhexa
c
c
c----------------------------------------------------------------------
c  TETRAHEDRON                                                        |
c----------------------------------------------------------------------
c
c  Type =
c
c     Linear : linear tetrahedron
c
c     TraTet : transfinite interpolation tetrahedron
c
c     CylTet : image of a linear tet through a global
c              system of cylindrical coordinates
c
c  VertNo - vertex points numbers
c
c  EdgeNo - edge curves numbers with a sign indicating orientation
c     +ncurve : consistency b/w curve (global) orientation and edge
c               local orientation
c     -ncurve : inconsistency
c
c  FigNo  - face figures nicknames; the orientations is the figure
c           (global) orientation wrt to the face local orientation
c     nick = ntrian*10 + norient [0,1,2,3,4,5]
c
c  Domain - domain number
c
c  Idata  - integer data for the tetrahedron
c
c----------------------------------------------------------------------
      type ttetra
        character(len=10)              :: Type
        integer                        :: VertNo(4)
        integer                        :: EdgeNo(6)
        integer                        :: FigNo(4)
        integer                        :: Domain
        integer, dimension(:), pointer :: Idata
      endtype ttetra
c
c
c----------------------------------------------------------------------
c  PYRAMID                                                            |
c----------------------------------------------------------------------
c
c  Type =
c
c     Linear  : pyramid with straight edges
c
c     TIPyram : transfinite interpolation pyramid
c
c  VertNo - vertex points numbers
c
c  EdgeNo - edge curves numbers with a sign indicating orientation
c     +ncurve : consistency b/w curve (global) orientation and edge
c               local orientation
c     -ncurve : inconsistency
c
c  FigNo  - face figures nicknames; the orientations is the figure
c           (global) orientation wrt to the face local orientation
c     lateral faces : nick = ntrian*10 + norient [0,1,2,3,4,5]
c     bottom  face  : nick = nrecta*10 + norient [0,1,2,3,4,5,6,7]
c
c  Domain - domain number
c
c  Idata  - integer data for the pyramid
c
c----------------------------------------------------------------------
      type ppyramid
        character(len=10)              :: Type
        integer                        :: VertNo(5)
        integer                        :: EdgeNo(8)
        integer                        :: FigNo(5)
        integer                        :: Domain
        integer, dimension(:), pointer :: Idata
      endtype ppyramid
c
c
c----------------------------------------------------------------------
c  ARRAYS OF 0,1,2,3-D OBJECTS                                        |
c----------------------------------------------------------------------
      type(surface  ), allocatable :: SURFACES(:)
      type(point    ), allocatable :: POINTS(:)
      type(point    ), allocatable :: POINTS_FG(:)
      type(ccurve   ), allocatable :: CURVES(:)
      type(triangle ), allocatable :: TRIANGLES(:)
      type(rectangle), allocatable :: RECTANGLES(:)
      type(hhexa    ), allocatable :: HEXAS(:)
      type(ttetra   ), allocatable :: TETRAS(:)
      type(pprism   ), allocatable :: PRISMS(:)
      type(ppyramid ), allocatable :: PYRAMIDS(:)
c
C$$$      interface input_DEFAULT
C$$$c
C$$$      subroutine input_DEFAULT_default
C$$$      end subroutine input_DEFAULT_default
C$$$c
C$$$      subroutine input_DEFAULT_from_file(fp)
C$$$c
C$$$      character(len=*), intent(in) :: fp
C$$$      end subroutine input_DEFAULT_from_file
C$$$      end interface
c
C$$$      interface input_NETGEN
c
C$$$      subroutine input_NETGEN_default
C$$$      end subroutine input_NETGEN_default
C$$$c
C$$$      subroutine input_NETGEN_from_file(fp)
C$$$c
C$$$      character(len=*), intent(in) :: fp
C$$$      end subroutine input_NETGEN_from_file
C$$$      end interface
c
      contains
c
c  ...allocates GMP data structures
      subroutine alloc_GMP
      integer :: i
c
      allocate(SURFACES(MAXSU))
      do i=1,MAXSU
        SURFACES(i)%Type = 'None'
        nullify(SURFACES(i)%Idata)
        nullify(SURFACES(i)%Rdata)
      enddo
c
      allocate(POINTS(MAXNP))
      do i=1,MAXNP
        POINTS(i)%Type = 'None'
        POINTS(i)%NrCurv = 0
        nullify(POINTS(i)%CurvNo)
        nullify(POINTS(i)%Idata)
        nullify(POINTS(i)%Rdata)
      enddo
c
      allocate(POINTS_FG(MAXNP_FG))
      do i=1,MAXNP_FG
        POINTS_FG(i)%Type = 'None'
        POINTS_FG(i)%NrCurv = 0
        nullify(POINTS_FG(i)%CurvNo)
        nullify(POINTS_FG(i)%Idata)
        nullify(POINTS_FG(i)%Rdata)
      enddo
c
      allocate(CURVES(MAXNC))
      do i=1,MAXNC
        CURVES(i)%Type = 'None'
        CURVES(i)%EndPoNo = 0
        CURVES(i)%NrFig = 0
        nullify(CURVES(i)%FigNo)
        nullify(CURVES(i)%Idata)
        nullify(CURVES(i)%Rdata)
      enddo
c
      allocate(TRIANGLES(MAXTR))
      do i=1,MAXTR
        TRIANGLES(i)%Type = 'None'
        TRIANGLES(i)%VertNo = 0
        TRIANGLES(i)%EdgeNo = 0
        TRIANGLES(i)%BlockNo = 0
        nullify(TRIANGLES(i)%Idata)
      enddo
c
      allocate(RECTANGLES(MAXRE))
      do i=1,MAXRE
        RECTANGLES(i)%Type = 'None'
        RECTANGLES(i)%VertNo = 0
        RECTANGLES(i)%EdgeNo = 0
        RECTANGLES(i)%BlockNo = 0
        nullify(RECTANGLES(i)%Idata)
        nullify(RECTANGLES(i)%Rdata)
      enddo
c
      allocate(PRISMS(MAXBT))
      do i=1,MAXBT
        PRISMS(i)%Type = 'None'
        PRISMS(i)%VertNo = 0
        PRISMS(i)%EdgeNo = 0
        PRISMS(i)%FigNo = 0
        PRISMS(i)%Domain = 0
      enddo
c
      allocate(HEXAS(MAXHE))
      do i=1,MAXHE
        HEXAS(i)%Type = 'None'
        HEXAS(i)%VertNo = 0
        HEXAS(i)%EdgeNo = 0
        HEXAS(i)%FigNo = 0
        HEXAS(i)%Domain = 0
        nullify(HEXAS(i)%Idata)
      enddo
c
      allocate(TETRAS(MAXTE))
      do i=1,MAXTE
        TETRAS(i)%Type = 'None'
        TETRAS(i)%VertNo = 0
        TETRAS(i)%EdgeNo = 0
        TETRAS(i)%FigNo = 0
        TETRAS(i)%Domain = 0
        nullify(TETRAS(i)%Idata)
      enddo
c
      allocate(PYRAMIDS(MAXPY))
      do i=1,MAXPY
        PYRAMIDS(i)%Type = 'None'
        PYRAMIDS(i)%VertNo = 0
        PYRAMIDS(i)%EdgeNo = 0
        PYRAMIDS(i)%FigNo = 0
        PYRAMIDS(i)%Domain = 0
        nullify(PYRAMIDS(i)%Idata)
      enddo
c
      end subroutine alloc_GMP
c
c
c----------------------------------------------------------------------
c
c  ...dump out GMP data structure
      subroutine dumpout_GMP
      integer :: iprint
      iprint=0
c
      ndump=31
      open(unit=ndump,file='files/dumpGMP',
     .     form='formatted',access='sequential',status='unknown')
c
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping GENERAL PARAMETRS'
      write(ndump,*)
     .           NDIM,MANDIM,
     .           NRSURFS,
     .           NRPOINT,NRCURVE,NRTRIAN,NRRECTA,
     .           NRHEXAS,NRTETRA,NRPRISM,NRPYRAM,
     .           NRDOMAIN
c
      write(ndump,*)
     .           MAXSU,
     .           MAXNP,MAXNC,MAXTR,MAXRE,
     .           MAXBT,MAXHE,MAXTE,MAXPY

c
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping SURFACES'
      do ns=1,MAXSU
        if (iprint.eq.1) write(*,*) 'dumpGMP: ns = ',ns
        write(ndump,*) SURFACES(ns)%Type
        if (iprint.eq.1) write(*,*) 'dumpGMP: SURFACES(ns)%Type = ',
     .                                        SURFACES(ns)%Type
c  .....Idata
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
c  .....Rdata
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
c
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
c  .....Idata
        if (associated(POINTS(np)%Idata)) then
          nn = ubound(POINTS(np)%Idata,1)
          nl = lbound(POINTS(np)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) POINTS(np)%Idata
        else
          write(ndump,*) 0
        endif
c  .....Rdata
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
c
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
c  .....Idata
        if (associated(CURVES(nc)%Idata)) then
          nn = ubound(CURVES(nc)%Idata,1)
          nl = lbound(CURVES(nc)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) CURVES(nc)%Idata
        else
          write(ndump,*) 0
        endif
c  .....Rdata
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
c
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping TRIANGLES'
      do nt=1,MAXTR
        write(ndump,*) TRIANGLES(nt)%Type
        write(ndump,*) TRIANGLES(nt)%VertNo(1:3)
        write(ndump,*) TRIANGLES(nt)%EdgeNo(1:3)
        write(ndump,*) TRIANGLES(nt)%BlockNo(1:2)
c  .....Idata
        if (associated(TRIANGLES(nt)%Idata)) then
          nn = ubound(TRIANGLES(nt)%Idata,1)
          nl = lbound(TRIANGLES(nt)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) TRIANGLES(nt)%Idata
        else
          write(ndump,*) 0
        endif
c  .....Rdata
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
c
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping RECTANGLES'
      do nr=1,MAXRE
        write(ndump,*) RECTANGLES(nr)%Type
        write(ndump,*) RECTANGLES(nr)%VertNo(1:4)
        write(ndump,*) RECTANGLES(nr)%EdgeNo(1:4)
        write(ndump,*) RECTANGLES(nr)%BlockNo(1:2)
c  .....Idata
        if (associated(RECTANGLES(nr)%Idata)) then
          nn = ubound(RECTANGLES(nr)%Idata,1)
          nl = lbound(RECTANGLES(nr)%Idata,1)
          write(ndump,*) nn
          write(ndump,*) nl
          write(ndump,*) RECTANGLES(nr)%Idata
        else
          write(ndump,*) 0
        endif
c  .....Rdata
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
c
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping PRISMS'
      do npri=1,MAXBT
        write(ndump,*) PRISMS(npri)%Type
        write(ndump,*) PRISMS(npri)%VertNo(1:6)
        write(ndump,*) PRISMS(npri)%EdgeNo(1:9)
        write(ndump,*) PRISMS(npri)%FigNo(1:5)
        write(ndump,*) PRISMS(npri)%Domain
c  .....Idata
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
c
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping HEXAS'
      do nh=1,MAXHE
        write(ndump,*) HEXAS(nh)%Type
        write(ndump,*) HEXAS(nh)%VertNo(1:8)
        write(ndump,*) HEXAS(nh)%EdgeNo(1:12)
        write(ndump,*) HEXAS(nh)%FigNo(1:6)
        write(ndump,*) HEXAS(nh)%Domain
c  .....Idata
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
c
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping TETRAS'
      do ntet=1,MAXTE
        write(ndump,*) TETRAS(ntet)%Type
        write(ndump,*) TETRAS(ntet)%VertNo(1:4)
        write(ndump,*) TETRAS(ntet)%EdgeNo(1:6)
        write(ndump,*) TETRAS(ntet)%FigNo(1:4)
        write(ndump,*) TETRAS(ntet)%Domain
c  .....Idata
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
c
      if (iprint.eq.1) write(*,*) 'dumpGMP: dumping PYRAMIDS'
      do npyr=1,MAXPY
        write(ndump,*) PYRAMIDS(npyr)%Type
        write(ndump,*) PYRAMIDS(npyr)%VertNo(1:5)
        write(ndump,*) PYRAMIDS(npyr)%EdgeNo(1:8)
        write(ndump,*) PYRAMIDS(npyr)%FigNo(1:5)
        write(ndump,*) PYRAMIDS(npyr)%Domain
c  .....Idata
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
c
      if (iprint.eq.1) write(*,*) 'dumpGMP: closing file'
      close(ndump)
c
      end subroutine dumpout_GMP
c
c
c----------------------------------------------------------------------
c
c  ...dump in GMP datastructure in LEGACY format (before Feb '11)
c
      subroutine dumpin_GMP_LEGACY(Fp)
c
      character(len=*) :: Fp
      iprint=0
c
      ndump=31
      open(unit=ndump,file=Fp,
     .     form='formatted',access='sequential',status='unknown')
c
      read(ndump,*)
     .           NDIM,MANDIM,
     .           NRSURFS,
     .           NRPOINT,NRCURVE,NRTRIAN,NRRECTA,
     .           NRHEXAS,NRTETRA,NRPRISM,NRPYRAM,
     .           NRDOMAIN
c
      read(ndump,*)
     .           MAXSU,
     .           MAXNP,MAXNC,MAXTR,MAXRE,
     .           MAXBT,MAXHE,MAXTE,MAXPY
c
ccc      if ((MAXSU.ne.MAXSU_old).or.
ccc     .    (MAXNP.ne.MAXNP_old).or.
ccc     .    (MAXNC.ne.MAXNC_old).or.
ccc     .    (MAXTR.ne.MAXTR_old).or.
ccc     .    (MAXRE.ne.MAXRE_old).or.
ccc     .    (MAXBT.ne.MAXBT_old).or.
ccc     .    (MAXHE.ne.MAXHE_old).or.
ccc     .    (MAXTE.ne.MAXTE_old).or.
ccc     .    (MAXPY.ne.MAXPY_old)) then
ccc        write(*,*) 'dumpin_GMP: INCOMPATIBLE GMP DATA STRUCTURES'
ccc        write(*,*)'MAXSU_old =',MAXSU_old
ccc        write(*,*)'MAXNP_old =',MAXNP_old
ccc        write(*,*)'MAXNC_old =',MAXNC_old
ccc        write(*,*)'MAXTR_old =',MAXTR_old
ccc        write(*,*)'MAXRE_old =',MAXRE_old
ccc        write(*,*)'MAXBT_old =',MAXBT_old
ccc        write(*,*)'MAXHE_old =',MAXHE_old
ccc        write(*,*)'MAXTE_old =',MAXTE_old
ccc        write(*,*)'MAXPY_old =',MAXPY_old
ccc        stop 1
ccc      endif
c
      if (iprint.eq.1) then
        write(*,*)'---------------------------------------------'
        write(*,*)'dumpin_GMP_LEGACY: '
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
      call alloc_GMP
c
      do ns=1,MAXSU
        read(ndump,*) SURFACES(ns)%Type
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(SURFACES(ns)%Idata(nn))
          read(ndump,*) SURFACES(ns)%Idata
        else
          nullify(SURFACES(ns)%Idata)
        endif
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(SURFACES(ns)%Rdata(nn))
          read(ndump,*) SURFACES(ns)%Rdata
        else
          nullify(SURFACES(ns)%Rdata)
        endif
      enddo
c
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
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(POINTS(np)%Idata(nn))
          read(ndump,*) POINTS(np)%Idata
        else
          nullify(POINTS(np)%Idata)
        endif
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(POINTS(np)%Rdata(nn))
          read(ndump,*) POINTS(np)%Rdata
        else
          nullify(POINTS(np)%Rdata)
        endif
      enddo
c
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
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(CURVES(nc)%Idata(nn))
          read(ndump,*) CURVES(nc)%Idata
        else
          nullify(CURVES(nc)%Idata)
        endif
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(CURVES(nc)%Rdata(nn))
          read(ndump,*) CURVES(nc)%Rdata
        else
          nullify(CURVES(nc)%Rdata)
        endif
      enddo
c
      do nt=1,MAXTR
        read(ndump,*) TRIANGLES(nt)%Type
        read(ndump,*) TRIANGLES(nt)%VertNo(1:3)
        read(ndump,*) TRIANGLES(nt)%EdgeNo(1:3)
        read(ndump,*) TRIANGLES(nt)%BlockNo(1:2)
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(TRIANGLES(nt)%Idata(nn))
          read(ndump,*) TRIANGLES(nt)%Idata
        else
          nullify(TRIANGLES(nt)%Idata)
        endif
      enddo
c
      do nr=1,MAXRE
        read(ndump,*) RECTANGLES(nr)%Type
        read(ndump,*) RECTANGLES(nr)%VertNo(1:4)
        read(ndump,*) RECTANGLES(nr)%EdgeNo(1:4)
        read(ndump,*) RECTANGLES(nr)%BlockNo(1:2)
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(RECTANGLES(nr)%Idata(nn))
          read(ndump,*) RECTANGLES(nr)%Idata
        else
          nullify(RECTANGLES(nr)%Idata)
        endif
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(RECTANGLES(nr)%Rdata(nn))
          read(ndump,*) RECTANGLES(nr)%Rdata
        else
          nullify(RECTANGLES(nr)%Rdata)
        endif
      enddo
c
      do npri=1,MAXBT
        read(ndump,*) PRISMS(npri)%Type
        read(ndump,*) PRISMS(npri)%VertNo(1:6)
        read(ndump,*) PRISMS(npri)%EdgeNo(1:9)
        read(ndump,*) PRISMS(npri)%FigNo(1:5)
        read(ndump,*) PRISMS(npri)%Domain
      enddo
c
      do nh=1,MAXHE
        read(ndump,*) HEXAS(nh)%Type
        read(ndump,*) HEXAS(nh)%VertNo(1:8)
        read(ndump,*) HEXAS(nh)%EdgeNo(1:12)
        read(ndump,*) HEXAS(nh)%FigNo(1:6)
        read(ndump,*) HEXAS(nh)%Domain
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(HEXAS(nh)%Idata(nn))
          read(ndump,*) HEXAS(nh)%Idata
        else
          nullify(HEXAS(nh)%Idata)
        endif
      enddo
c
      do ntet=1,MAXTE
        read(ndump,*) TETRAS(ntet)%Type
        read(ndump,*) TETRAS(ntet)%VertNo(1:4)
        read(ndump,*) TETRAS(ntet)%EdgeNo(1:6)
        read(ndump,*) TETRAS(ntet)%FigNo(1:4)
        read(ndump,*) TETRAS(ntet)%Domain
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(TETRAS(ntet)%Idata(nn))
          read(ndump,*) TETRAS(ntet)%Idata
        else
          nullify(TETRAS(ntet)%Idata)
        endif
      enddo
c
      do npyr=1,MAXPY
        read(ndump,*) PYRAMIDS(npyr)%Type
        read(ndump,*) PYRAMIDS(npyr)%VertNo(1:5)
        read(ndump,*) PYRAMIDS(npyr)%EdgeNo(1:8)
        read(ndump,*) PYRAMIDS(npyr)%FigNo(1:5)
        read(ndump,*) PYRAMIDS(npyr)%Domain
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(PYRAMIDS(npyr)%Idata(nn))
          read(ndump,*) PYRAMIDS(npyr)%Idata
        else
          nullify(PYRAMIDS(npyr)%Idata)
        endif
      enddo
c
      close(ndump)
c
      end subroutine dumpin_GMP_LEGACY
c
c
c----------------------------------------------------------------------
c
c  ...dump in GMP data structure
      subroutine dumpin_GMP(Fp)
c
      character(len=*), intent(in)  :: Fp
      integer, parameter :: ndump = 31
      iprint=0
c
      open(unit=ndump,file=Fp,
     .     form='formatted',access='sequential',status='unknown')
c
      read(ndump,*)
     .           NDIM,MANDIM,
     .           NRSURFS,
     .           NRPOINT,NRCURVE,NRTRIAN,NRRECTA,
     .           NRHEXAS,NRTETRA,NRPRISM,NRPYRAM,
     .           NRDOMAIN
c
      read(ndump,*)
     .           MAXSU,
     .           MAXNP,MAXNC,MAXTR,MAXRE,
     .           MAXBT,MAXHE,MAXTE,MAXPY
c
ccc      if ((MAXSU.ne.MAXSU_old).or.
ccc     .    (MAXNP.ne.MAXNP_old).or.
ccc     .    (MAXNC.ne.MAXNC_old).or.
ccc     .    (MAXTR.ne.MAXTR_old).or.
ccc     .    (MAXRE.ne.MAXRE_old).or.
ccc     .    (MAXBT.ne.MAXBT_old).or.
ccc     .    (MAXHE.ne.MAXHE_old).or.
ccc     .    (MAXTE.ne.MAXTE_old).or.
ccc     .    (MAXPY.ne.MAXPY_old)) then
ccc        write(*,*) 'dumpin_GMP: INCOMPATIBLE GMP DATA STRUCTURES'
ccc        write(*,*)'MAXSU_old =',MAXSU_old
ccc        write(*,*)'MAXNP_old =',MAXNP_old
ccc        write(*,*)'MAXNC_old =',MAXNC_old
ccc        write(*,*)'MAXTR_old =',MAXTR_old
ccc        write(*,*)'MAXRE_old =',MAXRE_old
ccc        write(*,*)'MAXBT_old =',MAXBT_old
ccc        write(*,*)'MAXHE_old =',MAXHE_old
ccc        write(*,*)'MAXTE_old =',MAXTE_old
ccc        write(*,*)'MAXPY_old =',MAXPY_old
ccc        stop 1
ccc      endif
c
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
c
      call alloc_GMP
c
      do ns=1,MAXSU
        read(ndump,*) SURFACES(ns)%Type
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(SURFACES(ns)%Idata(nl:nn))
          read(ndump,*) SURFACES(ns)%Idata
        else
          nullify(SURFACES(ns)%Idata)
        endif
c  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(SURFACES(ns)%Rdata(nl:nn))
          read(ndump,*) SURFACES(ns)%Rdata
        else
          nullify(SURFACES(ns)%Rdata)
        endif
      enddo
c
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
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(POINTS(np)%Idata(nl:nn))
          read(ndump,*) POINTS(np)%Idata
        else
          nullify(POINTS(np)%Idata)
        endif
c  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(POINTS(np)%Rdata(nl:nn))
          read(ndump,*) POINTS(np)%Rdata
        else
          nullify(POINTS(np)%Rdata)
        endif
      enddo
c
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
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(CURVES(nc)%Idata(nl:nn))
          read(ndump,*) CURVES(nc)%Idata
        else
          nullify(CURVES(nc)%Idata)
        endif
c  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(CURVES(nc)%Rdata(nl:nn))
          read(ndump,*) CURVES(nc)%Rdata
        else
          nullify(CURVES(nc)%Rdata)
        endif
      enddo
c
      do nt=1,MAXTR
        read(ndump,*) TRIANGLES(nt)%Type
        read(ndump,*) TRIANGLES(nt)%VertNo(1:3)
        read(ndump,*) TRIANGLES(nt)%EdgeNo(1:3)
        read(ndump,*) TRIANGLES(nt)%BlockNo(1:2)
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(TRIANGLES(nt)%Idata(nl:nn))
          read(ndump,*) TRIANGLES(nt)%Idata
        else
          nullify(TRIANGLES(nt)%Idata)
        endif
c  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(TRIANGLES(nt)%Rdata(nl:nn))
          read(ndump,*) TRIANGLES(nt)%Rdata
        else
          nullify(TRIANGLES(nt)%Rdata)
        endif
      enddo
c
      do nr=1,MAXRE
        read(ndump,*) RECTANGLES(nr)%Type
        read(ndump,*) RECTANGLES(nr)%VertNo(1:4)
        read(ndump,*) RECTANGLES(nr)%EdgeNo(1:4)
        read(ndump,*) RECTANGLES(nr)%BlockNo(1:2)
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(RECTANGLES(nr)%Idata(nl:nn))
          read(ndump,*) RECTANGLES(nr)%Idata
        else
          nullify(RECTANGLES(nr)%Idata)
        endif
c  .....Rdata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(RECTANGLES(nr)%Rdata(nl:nn))
          read(ndump,*) RECTANGLES(nr)%Rdata
        else
          nullify(RECTANGLES(nr)%Rdata)
        endif
      enddo
c
      do npri=1,MAXBT
        read(ndump,*) PRISMS(npri)%Type
        read(ndump,*) PRISMS(npri)%VertNo(1:6)
        read(ndump,*) PRISMS(npri)%EdgeNo(1:9)
        read(ndump,*) PRISMS(npri)%FigNo(1:5)
        read(ndump,*) PRISMS(npri)%Domain
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nn
          allocate(PRISMS(npri)%Idata(nl:nn))
          read(ndump,*) PRISMS(npri)%Idata
        else
          nullify(PRISMS(npri)%Idata)
        endif
      enddo
c
      do nh=1,MAXHE
        read(ndump,*) HEXAS(nh)%Type
        read(ndump,*) HEXAS(nh)%VertNo(1:8)
        read(ndump,*) HEXAS(nh)%EdgeNo(1:12)
        read(ndump,*) HEXAS(nh)%FigNo(1:6)
        read(ndump,*) HEXAS(nh)%Domain
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(HEXAS(nh)%Idata(nl:nn))
          read(ndump,*) HEXAS(nh)%Idata
        else
          nullify(HEXAS(nh)%Idata)
        endif
      enddo
c
      do ntet=1,MAXTE
        read(ndump,*) TETRAS(ntet)%Type
        read(ndump,*) TETRAS(ntet)%VertNo(1:4)
        read(ndump,*) TETRAS(ntet)%EdgeNo(1:6)
        read(ndump,*) TETRAS(ntet)%FigNo(1:4)
        read(ndump,*) TETRAS(ntet)%Domain
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(TETRAS(ntet)%Idata(nl:nn))
          read(ndump,*) TETRAS(ntet)%Idata
        else
          nullify(TETRAS(ntet)%Idata)
        endif
      enddo
c
      do npyr=1,MAXPY
        read(ndump,*) PYRAMIDS(npyr)%Type
        read(ndump,*) PYRAMIDS(npyr)%VertNo(1:5)
        read(ndump,*) PYRAMIDS(npyr)%EdgeNo(1:8)
        read(ndump,*) PYRAMIDS(npyr)%FigNo(1:5)
        read(ndump,*) PYRAMIDS(npyr)%Domain
c  .....Idata
        read(ndump,*) nn
        if (nn.gt.0) then
          read(ndump,*) nl
          allocate(PYRAMIDS(npyr)%Idata(nl:nn))
          read(ndump,*) PYRAMIDS(npyr)%Idata
        else
          nullify(PYRAMIDS(npyr)%Idata)
        endif
      enddo
c
      close(ndump)
c
      end subroutine dumpin_GMP
c
      subroutine print_GMP
c
      use control
#include "syscom.blk"
c
 10   continue
      write(*,*) 'print_GMP: SELECT:'
      write(*,*) 'EXIT....................................0'
      write(*,*) 'GENERAL DATA STRUCTURE PARAMETRS........1'
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
c
      select case(idec)
c
      case(0)
        return
c
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
        write(*,7006) MAXSU,
     .                MAXNP,MAXNC,MAXTR,MAXRE,
     .                MAXBT,MAXHE,MAXTE,MAXPY
 7006   format(' MAXSU,MAXNP,MAXNC,MAXTR,MAXRE,MAXBT,MAXHE,MAXTE,'
     .         'MAXPY = ',9i7)
c
      case(2)
        write(*,*) 'SET POINT NUMBER'
        read(*,*) np
        if ((np.le.0).or.(np.gt.NRPOINT)) go to 10
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
c
      case(3)
        write(*,*) 'SET CURVE NUMBER'
        read(*,*) nc
        if ((nc.le.0).or.(nc.gt.NRCURVE)) go to 10
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
c
      case(4)
        write(*,*) 'SET TRIANGLE NUMBER'
        read(*,*) nt
        if ((nt.le.0).or.(nt.gt.NRTRIAN)) go to 10
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
c
      case(5)
        write(*,*) 'SET RECTANGLE NUMBER'
        read(*,*) nr
        if ((nr.le.0).or.(nr.gt.NRRECTA)) go to 10
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
c
      case(6)
        write(*,*) 'SET PRISM NUMBER'
        read(*,*) npri
        if ((npri.le.0).or.(npri.gt.NRPRISM)) go to 10
        write(*,7701) npri,PRISMS(npri)%Type
 7701   format(' PRISM = ',i8,' TYPE = ',a10)
        write(*,7602) PRISMS(npri)%VertNo(1:6)
        write(*,7603) PRISMS(npri)%EdgeNo(1:9)
        write(*,7604) PRISMS(npri)%FigNo(1:5)
        write(*,7605) PRISMS(npri)%Domain
c
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
c
      case(8)
        write(*,*) 'SET TETRAHEDRON NUMBER'
        read(*,*) ntet
        if ((ntet.le.0).or.(ntet.gt.NRTETRA)) go to 10
        write(*,7801) ntet,TETRAS(ntet)%Type
 7801   format(' TETRA = ',i8,' TYPE = ',a10)
        write(*,7602) TETRAS(ntet)%VertNo(1:4)
        write(*,7603) TETRAS(ntet)%EdgeNo(1:6)
        write(*,7604) TETRAS(ntet)%FigNo(1:4)
        write(*,7605) TETRAS(ntet)%Domain
        if (associated(TETRAS(ntet)%Idata)) then
          write(*,8002) TETRAS(ntet)%Idata
        endif
c
      case(9)
        write(*,*) 'SET PYRAMID NUMBER'
        read(*,*) npyr
        if ((npyr.le.0).or.(npyr.gt.NRPYRAM)) go to 10
        write(*,7901) npyr,PYRAMIDS(npyr)%Type
 7901   format(' PYRAMID = ',i8,' TYPE = ',a10)
        write(*,7602) PYRAMIDS(npyr)%VertNo(1:5)
        write(*,7603) PYRAMIDS(npyr)%EdgeNo(1:8)
        write(*,7604) PYRAMIDS(npyr)%FigNo(1:5)
        write(*,7605) PYRAMIDS(npyr)%Domain
        if (associated(PYRAMIDS(npyr)%Idata)) then
          write(*,8002) PYRAMIDS(npyr)%Idata
        endif
c
      case(10)
        write(*,*) 'SET SURFACES NUMBER'
        read(*,*) ns
        if ((ns.le.0).or.(ns.gt.NRSURFS)) go to 10
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
c
      case default
        goto 10
      end select
      goto 10
c
      end subroutine print_GMP
c
c
      subroutine set_gmp_parameters(NDIM_loc,
     .                              MANDIM_loc,
     .                              MAXSU_loc,
     .                              MAXNP_loc,
     .                              MAXNC_loc,
     .                              MAXTR_loc,
     .                              MAXRE_loc,
     .                              MAXBT_loc,
     .                              MAXHE_loc,
     .                              MAXTE_loc,
     .                              MAXPY_loc)
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
c
      integer, parameter :: iprint=0
c
c  ...GMP control parameters...........................................
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
c
c  ...fine grid parameters.............................................
      NRPO_FG  = 0
      MAXNP_FG = MAXNP_loc
c
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
c
      end subroutine set_gmp_parameters
c
c
      subroutine print_GMP_parameters
      use environment , only : QUIET_MODE
c
      IF (.NOT. QUIET_MODE) then
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
ccc        write(*,8061)NRPO_FG,MAXNP_FG
ccc 8061   format(' NRPO_FG  = ',i12)
        write(*,*)''
      ENDIF
c
      end subroutine print_GMP_parameters
c
c
      integer function Ndomain_gmp(IGMP_block)
c
      integer IGMP_block, no, lab
      call decode(IGMP_block, no,lab)
      select case(lab)
      case(1)
        Ndomain_gmp = PRISMS(no)%Domain
      case(2)
        Ndomain_gmp = HEXAS(no)%Domain
      case(3)
        Ndomain_gmp = TETRAS(no)%Domain
      case(4)
        Ndomain_gmp = PYRAMIDS(no)%Domain
      case default
        write(*,*) 'Ndomain_gmp: lab = ',lab
        stop
      endselect
      end function
c
      end module GMP
