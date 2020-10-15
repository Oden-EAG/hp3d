!>  Purpose :
!!    module defines various parameters used to parametric elements
!!
!! @rev Dec 10
module parameters
  use error
  implicit none
  integer :: iprint_parameter
  !
  !  ...dimension of the problem
  integer, parameter :: NDIMEN = 3
  !
  integer, parameter :: NSTD_IN  = 5
  integer, parameter :: NSTD_OUT = 6
  integer, parameter :: NSTD_ERR = 0

  !  ...maximum number of physics attributes
  integer, parameter :: MAXPHYA=6
  !
  !  ...maximum order of approximation
  integer, parameter :: MAXP=4
  !
  !  ...modulo to encode polynomial orders
  !  ...(do not change unless using stand-alone shape functions package)
  integer, parameter :: MODORDER=10
  !
  !----------------------------------------------------------------------
  !
  !  ...maximum number of parent dof for a constrained dof
  integer, parameter :: NACDIM=(MAXP+1)**2
  !
  !  ...maximum number of nodes for a modified element
  integer, parameter :: MAXNODM=50
  !
  !  ...maximum number of generation
  integer, parameter :: MAXGEN = 20
  integer, parameter :: MAXQUEUE = 50
  !
  !  ...maximum number of quadrature points for 1D,2D,3D elements
  integer, parameter :: MAX_NINT1 =  MAXP+1
  integer, parameter :: MAX_NINT2 = (MAXP+1)**2
  integer, parameter :: MAX_NINT3 = (MAXP+1)**2*(MAXP+2)
  !
  !----------------------------------------------------------------------
  !  === ELEMENT ====
  !  ...max number of local dof for a 2D quad element
  integer, parameter :: MAXquadH = (MAXP+1)**2
  integer, parameter :: MAXquadE = 2*MAXP*(MAXP+1)
  integer, parameter :: MAXquadV = MAXquadE
  integer, parameter :: MAXquadQ = MAXP**2
  !
  !  ...max number of local dof for a 2D triangular element
  integer, parameter :: MAXtriaH = (MAXP+1)*(MAXP+2)/2
  integer, parameter :: MAXtriaE = MAXP*(MAXP+2)
  integer, parameter :: MAXtriaV = MAXtriaE
  integer, parameter :: MAXtriaQ = MAXP*(MAXP+1)/2
  !
  !----------------------------------------------------------------------
  !
  !  ...max number of local dof for a 3D brick element
  integer, parameter :: MAXbrickH = (MAXP+1)**3
  integer, parameter :: MAXbrickE = 3*MAXP*(MAXP+1)**2
  integer, parameter :: MAXbrickV = 3*MAXP**2*(MAXP+1)
  integer, parameter :: MAXbrickQ = MAXP**3
  !
  !  ...max number of local dof for a 3D prism element
  integer, parameter :: MAXprismH = MAXtriaH*(MAXP+1)
  integer, parameter :: MAXprismE = MAXtriaE*(MAXP+1) + MAXtriaH*MAXP
  integer, parameter :: MAXprismV = MAXtriaV*MAXP + MAXtriaQ*(MAXP+1)
  integer, parameter :: MAXprismQ = MAXtriaQ*MAXP
  !
  !  ...max number of local dof for a 3D tetrahedral  element
  integer, parameter :: MAXtetraH = (MAXP+1)*(MAXP+2)*(MAXP+3)/6
  integer, parameter :: MAXtetraE = MAXP*(MAXP+2)*(MAXP+3)/2
  integer, parameter :: MAXtetraV = MAXP*(MAXP+1)*(MAXP+3)/2
  integer, parameter :: MAXtetraQ = MAXP*(MAXP+1)*(MAXP+2)/6
  !
  !  ...max number of local dof for a 3D pyramid  element
  integer, parameter :: MAXpyramH = 5+8*(MAXP-1)+(MAXP-1)**2+(MAXP-2)*(MAXP-1)*2+(MAXP-1)**3
  integer, parameter :: MAXpyramE = 8*MAXP+2*MAXP*(MAXP-1)+4*MAXP*(MAXP-1)+3*(MAXP-1)**2*MAXP
  integer, parameter :: MAXpyramV = MAXquadQ+4*MAXtriaQ+3*(MAXP-1)*MAXP**2
  integer, parameter :: MAXpyramQ = MAXP**3
  !
  !----------------------------------------------------------------------
  !  ==== NODE ===
  !  ...maximum number of dof for 'mdlt'
  integer, parameter :: MAXmdltH=(MAXP-2)*(MAXP-1)/2
  integer, parameter :: MAXmdltE=(MAXP-1)*MAXP
  integer, parameter :: MAXmdltV=MAXmdltE
  integer, parameter :: MAXmdltQ=MAXP*(MAXP+1)/2
  !
  !  ...maximum number of dof for 'mdlq'
  integer, parameter :: MAXmdlqH=(MAXP-1)**2
  integer, parameter :: MAXmdlqE=2*MAXP*(MAXP-1)
  integer, parameter :: MAXmdlqV=MAXmdlqE
  integer, parameter :: MAXmdlqQ=MAXP**2
  !
  !  ...maximum number of dof for 'mdlb'
  integer, parameter :: MAXmdlbH=(MAXP-1)**3
  integer, parameter :: MAXmdlbE=3*MAXP*(MAXP-1)**2
  integer, parameter :: MAXmdlbV=3*MAXP**2*(MAXP-1)
  integer, parameter :: MAXmdlbQ=MAXbrickQ
  !
  !  ...maximum number of dof for 'mdln'
  integer, parameter :: MAXmdlnH=(MAXP-3)*(MAXP-2)*(MAXP-1)/6
  integer, parameter :: MAXmdlnE=(MAXP-2)*(MAXP-1)*MAXP/2
  integer, parameter :: MAXmdlnV=(MAXP-1)*MAXP*(MAXP+1)/2
  integer, parameter :: MAXmdlnQ=MAXtetraQ
  !
  !  ...maximum number of dof for 'mdlp'
  integer, parameter :: MAXmdlpH=MAXmdltH*(MAXP-1)
  integer, parameter :: MAXmdlpE=MAXmdltE*(MAXP-1)+MAXmdltH*MAXP
  integer, parameter :: MAXmdlpV=MAXmdltV*MAXP+MAXmdltQ*(MAXP-1)
  integer, parameter :: MAXmdlpQ=MAXprismQ
  !
  !  ...maximum number of dof for 'mdld'
  integer, parameter :: MAXmdldH=(MAXP-1)**3
  integer, parameter :: MAXmdldE=3*(MAXP-1)**2*MAXP
  integer, parameter :: MAXmdldV=3*(MAXP-1)*MAXP**2
  integer, parameter :: MAXmdldQ=MAXpyramQ
  !
  !----------------------------------------------------------------------

  !  ...number of maximum timesteps
  ! integer, save :: MAXTIMESTEP

  !  ...number of copies of variables stored
  integer, save :: NRCOMS
  !
  !  ...maximum number of right-hand sides (load vectors)
  integer, save :: MAXNRHS
  !
  !  ...maximum number of solution H1 components for all possible cases
  integer, save :: MAXEQNH
  !
  !  ...maximum number of solution components discretized
  !     with H(curl)-conforming shape functions
  integer, save :: MAXEQNE
  !
  !  ...maximum number of solution components discretized
  !     with H(div)-conforming shape functions
  integer, save :: MAXEQNV
  !
  !  ...maximum number of solution components discretized
  !     with L^2-conforming shape functions
  integer, save :: MAXEQNQ
  !
  !----------------------------------------------------------------------
  ! real parameters
#if C_MODE
  complex(8), parameter :: ZERO = (0.d0,0.d0)
  complex(8), parameter :: ZONE = (1.d0,0.d0)
  complex(8), parameter :: ZEYE = (1.d0,0.d0)
  complex(8), parameter :: ZIMG = (0.d0,1.d0)
#else
  real(8)   , parameter :: ZERO = 0.d0
  real(8)   , parameter :: ZONE = 1.d0
  complex(8), parameter :: ZIMG = (0.d0,1.d0)
#endif
  !
  contains

    subroutine set_parameters(LNrcoms, LMaxnrhs, &
         LMaxeqnh, LMaxeqne, LMaxeqnv, LMaxeqnq)
      integer :: LNrcoms, LMaxnrhs, &
           LMaxeqnh, LMaxeqne, LMaxeqnv, LMaxeqnq
      NRCOMS  = LNrcoms
      MAXNRHS = LMaxnrhs

      MAXEQNH = LMaxeqnh;      MAXEQNE = LMaxeqne
      MAXEQNV = LMaxeqnv;      MAXEQNQ = LMaxeqnq
    end subroutine set_parameters
end module parameters
