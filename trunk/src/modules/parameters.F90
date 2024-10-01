!> @brief Defines various parameters for parametric elements
!> @date  Sep 2023
module parameters
!
   use error
!
   implicit none
!
!..dimension of the problem
   integer, parameter :: NDIMEN = 3
!
   integer, parameter :: NSTD_IN  = 5
   integer, parameter :: NSTD_OUT = 6
   integer, parameter :: NSTD_ERR = 0
!
!..maximum order of approximation: 1,...,9
   integer, parameter :: MAXP = 6
!
!..modulo to encode polynomial orders
!  (do not change unless using stand-alone shape functions package)
   integer, parameter :: MODORDER = 10
!
!----------------------------------------------------------------------
!
!..maximum number of parent dof for a constrained dof
   integer, parameter :: NACDIM = (MAXP+1)**2
!
!..maximum number of nodes for a modified element
   integer, parameter :: MAXNODM = 50
!
!..maximum number of generation
   integer, parameter :: MAXGEN = 20
   integer, parameter :: MAXQUEUE = 50
!
!..maximum number of quadrature points for 1D,2D,3D elements
   integer, parameter :: MAX_NINT1 =  MAXP+1
   integer, parameter :: MAX_NINT2 = (MAXP+1)**2
   integer, parameter :: MAX_NINT3 = (MAXP+1)**2*(MAXP+2)
!
!----------------------------------------------------------------------
!  === ELEMENT ====
!..max number of local dof for a 2D quad element
   integer, parameter :: MAXquadH = (MAXP+1)**2
   integer, parameter :: MAXquadE = 2*MAXP*(MAXP+1)
   integer, parameter :: MAXquadV = MAXquadE
   integer, parameter :: MAXquadQ = MAXP**2
!
!..max number of local dof for a 2D triangular element
   integer, parameter :: MAXtriaH = (MAXP+1)*(MAXP+2)/2
   integer, parameter :: MAXtriaE = MAXP*(MAXP+2)
   integer, parameter :: MAXtriaV = MAXtriaE
   integer, parameter :: MAXtriaQ = MAXP*(MAXP+1)/2
!
!----------------------------------------------------------------------
!
!..max number of local dof for a 3D brick element
   integer, parameter :: MAXbrickH = (MAXP+1)**3
   integer, parameter :: MAXbrickE = 3*MAXP*(MAXP+1)**2
   integer, parameter :: MAXbrickV = 3*MAXP**2*(MAXP+1)
   integer, parameter :: MAXbrickQ = MAXP**3
!
!..max number of local dof for a 3D prism element
   integer, parameter :: MAXprismH = MAXtriaH*(MAXP+1)
   integer, parameter :: MAXprismE = MAXtriaE*(MAXP+1) + MAXtriaH*MAXP
   integer, parameter :: MAXprismV = MAXtriaV*MAXP + MAXtriaQ*(MAXP+1)
   integer, parameter :: MAXprismQ = MAXtriaQ*MAXP
!
!..max number of local dof for a 3D tetrahedral element
   integer, parameter :: MAXtetraH = (MAXP+1)*(MAXP+2)*(MAXP+3)/6
   integer, parameter :: MAXtetraE = MAXP*(MAXP+2)*(MAXP+3)/2
   integer, parameter :: MAXtetraV = MAXP*(MAXP+1)*(MAXP+3)/2
   integer, parameter :: MAXtetraQ = MAXP*(MAXP+1)*(MAXP+2)/6
!
!..max number of local dof for a 3D pyramid element
   integer, parameter :: MAXpyramH = 5+8*(MAXP-1)+(MAXP-1)**2+(MAXP-2)*(MAXP-1)*2+(MAXP-1)**3
   integer, parameter :: MAXpyramE = 8*MAXP+2*MAXP*(MAXP-1)+4*MAXP*(MAXP-1)+3*(MAXP-1)**2*MAXP
   integer, parameter :: MAXpyramV = MAXquadQ+4*MAXtriaQ+3*(MAXP-1)*MAXP**2
   integer, parameter :: MAXpyramQ = MAXP**3
!
!----------------------------------------------------------------------
!  ==== NODE ===
!..maximum number of dof for 'mdlt' (tria middle node)
   integer, parameter :: MAXmdltH = (MAXP-2)*(MAXP-1)/2
   integer, parameter :: MAXmdltE = (MAXP-1)*MAXP
   integer, parameter :: MAXmdltV = MAXmdltE
   integer, parameter :: MAXmdltQ = MAXP*(MAXP+1)/2
!
!..maximum number of dof for 'mdlq' (quad middle node)
   integer, parameter :: MAXmdlqH = (MAXP-1)**2
   integer, parameter :: MAXmdlqE = 2*MAXP*(MAXP-1)
   integer, parameter :: MAXmdlqV = MAXmdlqE
   integer, parameter :: MAXmdlqQ = MAXP**2
!
!..maximum number of dof for 'mdlb' (bric middle node)
   integer, parameter :: MAXmdlbH = (MAXP-1)**3
   integer, parameter :: MAXmdlbE = 3*MAXP*(MAXP-1)**2
   integer, parameter :: MAXmdlbV = 3*MAXP**2*(MAXP-1)
   integer, parameter :: MAXmdlbQ = MAXbrickQ
!
!..maximum number of dof for 'mdln' (tetr middle node)
   integer, parameter :: MAXmdlnH = (MAXP-3)*(MAXP-2)*(MAXP-1)/6
   integer, parameter :: MAXmdlnE = (MAXP-2)*(MAXP-1)*MAXP/2
   integer, parameter :: MAXmdlnV = (MAXP-1)*MAXP*(MAXP+1)/2
   integer, parameter :: MAXmdlnQ = MAXtetraQ
!
!..maximum number of dof for 'mdlp' (pris middle node)
   integer, parameter :: MAXmdlpH = MAXmdltH*(MAXP-1)
   integer, parameter :: MAXmdlpE = MAXmdltE*(MAXP-1)+MAXmdltH*MAXP
   integer, parameter :: MAXmdlpV = MAXmdltV*MAXP+MAXmdltQ*(MAXP-1)
   integer, parameter :: MAXmdlpQ = MAXprismQ
!
!..maximum number of dof for 'mdld' (pyra middle node)
   integer, parameter :: MAXmdldH = (MAXP-1)**3
   integer, parameter :: MAXmdldE = 3*(MAXP-1)**2*MAXP
   integer, parameter :: MAXmdldV = 3*(MAXP-1)*MAXP**2
   integer, parameter :: MAXmdldQ = MAXpyramQ
!
!----------------------------------------------------------------------
!
!..number of solution component sets (COMS) stored
!  NRCOMS = 1 by default (only one copy of the solution DOFs)
   integer, save :: NRCOMS
!
!..solution component set being assembled and solved for
!  currently, only N_COMS=1 is supported
   integer, save :: N_COMS
!
!..number of right-hand sides (load vectors)
!  currently, only NRRHS=1 is supported
   integer, save :: NRRHS
!
!..maximum number of solution H1 components for all possible cases
!  MAXEQNH := NRHVAR * NRRHS
   integer, save :: MAXEQNH
!
!..maximum number of solution components discretized
!  with H(curl)-conforming shape functions
!  MAXEQNE := NREVAR * NRRHS
   integer, save :: MAXEQNE
!
!..maximum number of solution components discretized
!  with H(div)-conforming shape functions
!  MAXEQNV := NRVVAR * NRRHS
   integer, save :: MAXEQNV
!
!..maximum number of solution components discretized
!  with L^2-conforming shape functions
!  MAXEQNQ := NRQVAR * NRRHS
   integer, save :: MAXEQNQ
!
!..auxiliary parameter
   logical, private, save :: IS_SET_PARAMETERS = .false.
!
!----------------------------------------------------------------------
!..real parameters
   real(8)   , parameter :: rZERO = 0.d0
   real(8)   , parameter :: rONE  = 1.d0
#if HP3D_COMPLEX
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
!
contains
!
!----------------------------------------------------------------------
!> @brief   Set user-defined parameters
!> @date    Sep 2023
   subroutine set_parameters(LNrComs,LNrRhs)
!
      integer, intent(in) :: LNrComs,LNrRhs
!
      if (IS_SET_PARAMETERS) then
         write(*,1000) 'Parameters may only be set once.'
         stop
      endif
      if (LNrRhs .ne. 1) then
         write(*,1000) 'Currently, only NRRHS=1 is supported.'
         stop
      endif
      if (LNrComs .lt. 1) then
         write(*,2000) 'Invalid parameter: LNrComs',LNrComs
         stop
      endif
 1000 format("set_parameters: ",A," stop.")
 2000 format("set_parameters: ",A," = ",I9,"; stop.")
!
!  ...set number of solution copies and number of right-hand sides
      NRCOMS = LNrComs
      NRRHS  = LNrRhs
!
!  ...by default, the first solution component set is assembled and solved for
      N_COMS = 1
!
      IS_SET_PARAMETERS = .true.
!
   end subroutine set_parameters
!
end module parameters
