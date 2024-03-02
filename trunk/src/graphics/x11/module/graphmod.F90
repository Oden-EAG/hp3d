#if HP3D_USE_X11

!> @brief define the workspace for graphics
module graphmod
  !
  implicit none
  !
  ! ** coltab.blk :: color table
  integer, parameter :: NR_COLORS = 200
  integer, dimension(NR_COLORS) :: NPCOL, IRED, IBLUE, IGREEN
  !
  !    array of colors NPCOL:
  !        NPCOL(1)-white,NPCOL(2)-black,
  !        NPCOL(3-10)-low-order-to-high-order-p-approx
  !        NPCOL(11-NR_COLORS) - spectrum of colors
  !    arrays IRED,IBLUE,IGREEN contain mixing specifications
  !           for the defined colors
  !
  !    parameter NR_COLORS depends upon your graphics card ans it
  !    should be determined experimentally by calling routine
  !    graph_util/colorset with printing flag up and determining
  !    how many colors can be supported

  ! ** gwind.blk windows parameter
  integer :: ISIZE(4),IWINDNUM,IDISPLAY_TYPE
  ! modify IWINDL, IWINDH to change X11 window size:
  integer :: IWINDL = 1334, IWINDH = 800
  real(8) :: RWINDL,RWINDH,RMARGIN,XLENGTH,YLENGTH,RSIZE(4),RANGE(4)
  !      ISIZE   - window size
  !      ISIZE(1) = 0
  !      ISIZE(2) = 0
  !      ISIZE(3) = IWINDL
  !      ISIZE(4) = IWINDH
  !      IWINDNUM - window number
  !      RSIZE - different for different displays
  !      RWINDH,RWINDL - as above
  !      RANGE - to change scale
  !      RMARGIN = .05d0*RWINDH
  !      XLENGTH = RWINDL - 3.d0*rmargin
  !      YLENGTH = RWINDH - 2.d0*rmargin
  !
  !      IDISPLAY_TYPE  -  type of display
  real(8) :: RN(3), RMTR(3,3)
  !
  !     RN  - components of the normal unit vector for a projection
  !           plane
  !     RMTR - matrix of transformation from cartesian to
  !           observers system

  ! ** graphsp.blk
  ! parameters
  integer, parameter :: MAXNRINVBL  = 20000
  integer, parameter :: MAXNRCURVBL = 1000
  integer, parameter :: MAXNRDOMAIN = 20

  ! list of invisible blocks
  integer, dimension(MAXNRINVBL)  :: IGINV
  real(8), dimension(300)         :: RTRMP
  integer, dimension(MAXNRCURVBL) :: NLINBLOCKS
  integer, dimension(MAXNRDOMAIN) :: NDOMAIN
  !
  ! physical attribute
  integer :: NRPHY_DISP

  real(8), dimension(1:3,1:2) :: BOX_CUT
  integer :: IBOX_CUT

  ! ** Visible objects :: gparams.blk
  integer :: NRVISTR, NRINVBL, NRCURVBL
  !
  !       NRVISTR - number of visible triangles
  !       NRINVBL - number of invisible blocks
  !       NRCURVBL - number of curvilinear blocks

  ! ** Color boxes :: gbox.blk
  real(8), dimension(1:2,4,0:8) :: XY_BOX

  ! ** Select the display quantity :: gselect.blk
  integer, save :: ISELECT

  ! ** Define sclae of triangles :: gscale.blk
  integer :: NRSUB
  real(8) :: DX, DIMOB(3), XCENTR(3), XCIM(2), XY(2,3), &
             SIZE, XCWIN(2), DIMIM, XEX(6), CLPL(4)

  !     NRSUB - number of subdivisions for each edge
  !     DX    = 1/NRSUB
  !     XCENTR- coordinates of the center of the object ( (xmax+xmin)/2 )
  !             in the physical space
  !     DIMOB = (xmax - xmin)/2
  !     XCIM  - coordinates of the center of the displayed part of object
  !     DIMIM - rescaled dimensions of the object
  !     XEX   - extreme values of coordinates for the object
  !     CLPL  - coefficients defining clipping plane for cross
  !             sections ('1*x+'2*y+'3*z+'4=0)
  !
  !  The physical object is projected onto a plane perpendicular to
  !  the direction of observation (point of view), and rescaled
  !  according to the formula:
  !
  !     xy(j) = (xy(j)-XCIM(j))/DIMIM*WSIZE + XCWIN(j)
  !
  !     where WSIZE - size of the window
  !           XCWIN - coordinates of the centre of the window
  !
  ! ** Workspace
  !--------------------------------------------------------------
  integer, save :: MXIGTR=0,MXIGSTR=0,MXRGTRZ=0
  integer, allocatable, dimension(:) :: IGTRCU, IGTRNO, IGTR, IGSTR
  real(8), allocatable, dimension(:) :: RGTR, RGTRZ
  logical :: INITIALIZED = .false.
  !  explanation of variables:
  !     RGTR   - real storage place for points coordinates
  !              for graphics (in observer's system) for
  !              triangles vertices
  !     IGTRCU - information on triangle sides, color etc
  !     IGTRNO - information about nodes numbers or small triangles
  !     IGSTR  - information about colour of small triangles
  !     RGTRZ  - depth-coordinates of triangles' mid-points
  !     IGTR   - ordered list of triangles (back-to-front)
  !     IGINV  - invisible blocks
end module graphmod

#else

module graphmod
   implicit none
end module graphmod

#endif
