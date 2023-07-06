!> Purpose : module stores a number of control parameters
module control
   save
!
   integer            :: SHAPE_FLAG
   integer, parameter :: PEANO_SHAPE     = 1
   integer, parameter :: LEGENDRE_SHAPE  = 2
!
!..a flag indicating whether the exact solution is known
   integer            :: NEXACT
!
!..a symmetry flag
   integer            :: ISYM_FLAG
   integer, parameter :: SYMMETRIC    = 1
   integer, parameter :: UNSYMMETRIC  = 2
!
!..static condensation flag
   logical            :: ISTC_FLAG
!
!..adaptive integration flags (LEGACY)
   integer            :: INT_LOAD
   integer            :: INT_NEUMANN
!
!..extension flag for supporting construction of the extension operator (LEGACY)
   integer            :: NFLAG2GRID
!
!..input mode flag
   integer            :: INPUT_FILE
   integer, parameter :: DUMPIN_LEGACY_  = -1
   integer, parameter :: DUMPIN_         =  0
   integer, parameter :: LEGACY_         =  1
   integer, parameter :: DEFAULT_        =  2
   integer, parameter :: COMPACT_        =  3
   integer, parameter :: RECONSTRUCT_    =  4
   integer, parameter :: NETGEN_         =  5
   integer, parameter :: LAGRANGE_       =  6
!
!..geometry tolerance
   real(8)            :: GEOM_TOL
!
!..overintegration index
   integer            :: INTEGRATION
!$OMP THREADPRIVATE (INTEGRATION)
!
!..exact geometry elements flag
   integer            :: EXGEOM
!
   interface read_control
      module procedure read_control_from_default
      module procedure read_control_from_file
   end interface
!
contains
!
!--------------------------------------------------------------------------
!> Purpose : routine reads in global control parameters from 'file/control'
   subroutine read_control_from_default
      call read_control_from_file('files/control')
   end subroutine read_control_from_default
!
!--------------------------------------------------------------------------
!> Purpose : routine reads in global control parameters from file
   subroutine read_control_from_file(fp)
      use environment , only : QUIET_MODE
      use paraview    , only : PARAVIEW_GEOM
      use MPI
!
      implicit none
!
!  ...input argument
      character(len=*), intent(in) :: fp
!
!  ...file handler for control
      integer :: ncontrol
!
!  ...a file containing various control parameters
      ncontrol=132
      open(unit=ncontrol,file=fp, &
         form='formatted',access='sequential',status='old',action='read')
!
      if (.not. QUIET_MODE) write(*,*)'-- hp3d Control Parameters --'
!
!  ...read in control parameters
      read(ncontrol,*) SHAPE_FLAG
      select case(SHAPE_FLAG)
         case(PEANO_SHAPE)
            if (.not. QUIET_MODE) write(*,*)'PEANO SHAPE FUNCTIONS USED'
         case(LEGENDRE_SHAPE)
            if (.not. QUIET_MODE) write(*,*)'INT. LEGENDRE POLYNOMIALS USED'
         case default
            write(*,*)'invalid value for SHAPE_FLAG!'
            stop
      endselect
!
      read(ncontrol,*) NEXACT
      select case(NEXACT)
         case(0)
            if (.not. QUIET_MODE) write(*,*)'EXACT SOLUTION UNKNOWN'
         case(1,2)
            if (.not. QUIET_MODE) write(*,*)'EXACT SOLUTION KNOWN'
         case default
            write(*,*)'invalid value for NEXACT!'
            stop
      endselect
!
      read(ncontrol,*) ISYM_FLAG
      select case(ISYM_FLAG)
         case(UNSYMMETRIC)
            if (.not. QUIET_MODE) write(*,*)'UNSYMMETRIC PROBLEM'
         case(SYMMETRIC)
            if (.not. QUIET_MODE) write(*,*)'SYMMETRIC PROBLEM'
         case default
            write(*,*)'invalid value for ISYM_FLAG!'
            stop
      endselect
!
      ISTC_FLAG = .false.
!
      read(ncontrol,*) INT_LOAD
      select case(INT_LOAD)
      case(0)
         if (.not. QUIET_MODE) write(*,*)'GAUSSIAN QUADRATURE (VOLUME LOAD)'
      case(1)
         if (.not. QUIET_MODE) write(*,*)'ADAPTIVE QUADRATURE (VOLUME LOAD)'
      endselect
!
      read(ncontrol,*) INT_NEUMANN
      select case(INT_NEUMANN)
         case(0)
            if (.not. QUIET_MODE) write(*,*)'GAUSSIAN QUADRATURE (NEUMANN BC)'
         case(1)
            if (.not. QUIET_MODE) write(*,*)'ADAPTIVE QUADRATURE (NEUMANN BC)'
      endselect
!
      read(ncontrol,*) NFLAG2GRID
      select case(NFLAG2GRID)
         case(0)
            if (.not. QUIET_MODE) write(*,*)'EXTENSION OPERATOR NOT SUPPORTED'
         case(1)
            if (.not. QUIET_MODE) write(*,*)'EXTENSION OPERATOR SUPPORTED'
      endselect
!
      read(ncontrol,*) INPUT_FILE
      select case(INPUT_FILE)
         case(DUMPIN_LEGACY_)
            if (.not. QUIET_MODE) write(*,*)'DUMPIN LEGACY FORMAT'
         case(DUMPIN_)
            if (.not. QUIET_MODE) write(*,*)'DUMPIN FORMAT'
         case(LEGACY_)
            if (.not. QUIET_MODE) write(*,*)'LEGACY FORMAT'
         case(DEFAULT_)
            if (.not. QUIET_MODE) write(*,*)'DEFAULT FORMAT'
         case(COMPACT_)
            if (.not. QUIET_MODE) write(*,*)'COMPACT FORMAT'
         case(RECONSTRUCT_)
            if (.not. QUIET_MODE) write(*,*)'RECONSTRUCTION FORMAT'
         case(NETGEN_)
            if (.not. QUIET_MODE) write(*,*)'NETGEN FORMAT'
         case(LAGRANGE_)
            if (.not. QUIET_MODE) write(*,*)'LAGRANGE FORMAT'
         case default
            write(*,*)'invalid value for INPUT_FILE!'
            stop
      endselect
!
      read(ncontrol,*) GEOM_TOL
      if (.not. QUIET_MODE) write(*,8000)GEOM_TOL
 8000 format(' GEOMETRY TOLERANCE = ',es11.4)
!
      read(ncontrol,*) INTEGRATION
      if (.not. QUIET_MODE) write(*,8001)INTEGRATION
 8001 format(' INTEGRATION        = ',i1)
!
      read(ncontrol,*) EXGEOM
      select case(EXGEOM)
         case(0)
            if (.not. QUIET_MODE) write(*,*)'ISOPARAMETRIC ELEMENTS'
         case(1)
            if (.not. QUIET_MODE) write(*,*)'EXACT GEOMETRY ELEMENTS'
         case default
            write(*,*)'invalid value for EXGEOM!'
            stop
      endselect
      if (.not. QUIET_MODE) write(*,*)''
!
      close (ncontrol)
!
!  ...initialize Paraview ex./iso. geometry flag
      PARAVIEW_GEOM = EXGEOM
!
   end subroutine read_control_from_file
!
end module control
