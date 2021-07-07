!----------------------------------------------------------------------------------
!> Purpose : test multiple polynomial projections for each problem
!!
!> @date : Nov 14
!----------------------------------------------------------------------------------
!
subroutine test_1
!      
      use PROJ        , only : NPX,NPY,NPZ,NUMBER_OF_RHS,ITAG,ITEST,ISOLVER, &
                               NEDELEC,ICOMP
      use physics     , only : NR_PHYSA
!
      implicit none
      integer :: itag0,idec
      integer, dimension(NR_PHYSA) :: iphys
!      
!----------------------------------------------------------------------------------
!
!     -- INTERACTIVE MODE --
      if (ITEST == 0) then
!
!       need to compute error for ALL physical attributes
        iphys(1:NR_PHYSA)=1
!        
!       force loop entry
        idec=1
!
!       repeat in infinite loop
        do while (idec /= 0)
!        
!         give user chance to choose solver
          write(*,*)'Solve: 1 - Frontal ; 2 - MUMPS'
          read( *,*) ISOLVER
!
!         chose b/w single projections and all problem-dependent projections
          write(*,*)'Test: 0 - Single projection ; 1 - All projections'
          read( *,*) idec
!
!         -- SINGLE PROJECTION --
          if (idec == 0) then
  
!           set monomial
            write(*,*) 'Set: NPX, NPY, NPZ ='
            read( *,*) NPX,NPY,NPZ
!            
!           create tag
            itag0 = ITAG*10000 + NPX*100 + NPY*10 + NPZ*1
!       
!           force "exact" to return lower order polynomials for all spaces
            NEDELEC=0 ; ICOMP=0
!
!           solve
            select case(ISOLVER)
            case(1) ; call solve1(         NUMBER_OF_RHS)
            case(2) ; call mumps_sc('G')
!!            case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
            endselect
!
!           compute error
            call compute_error(iphys,itag0)
!
!         -- MULTIPLE PROJECTIONS --
          else ; call test_1_aux
          endif
!
          write(*,*)'Test 1: 0 - Return ; 1 - Continue testing'
          read( *,*) idec
!
        enddo
!
!     -- TESTING MODE
      else ; call test_1_aux
      endif
!
!
endsubroutine test_1
!
!
!
!----------------------------------------------------------------------------------
!
subroutine test_1_aux
!
      use PROJ        , only : NPX,NPY,NPZ,IORDER_TETR,IORDER_BRIC, & 
                                           IORDER_PRIS,IORDER_PYRA, &
                                           PROBLEM,NUMBER_OF_RHS  , &
                                           NEDELEC, ICOMP, ITAG, ISOLVER
      use environment , only : QUIET_MODE
      use physics     , only : NR_PHYSA
!!      use uhm2        , only : UHM_SOLVER_PTR
!
      implicit none
      integer :: i,iorder,ispace,p,ncomp,itag0,idec
      integer, dimension(NR_PHYSA) :: iphys
!
!----------------------------------------------------------------------------------
!
!     need to compute error for ALL physical attributes
      iphys(1:NR_PHYSA)=1
!
!     -- TESTING MODE --
      select case(PROBLEM)
!
!----------------------------------------------------------------------------------
!     TETRAHEDRON                                                                 |
!     HYBRID : PRISM + HEXA + TET (+ PYRAMID)                                     |
!----------------------------------------------------------------------------------
      case('tetr','hybr')
!
!       L O W E R    O R D E R    P O L Y N O M I A L S
!
!       force "exact" to return lower order polynomials
        NEDELEC=0 ; ICOMP=0
!        
        do NPX=0,IORDER_TETR-1
        do NPY=0,IORDER_TETR-1
        do NPZ=0,IORDER_TETR-1
!        
          if (NPX+NPY+NPZ > IORDER_TETR-1)  cycle
!  
!         solve
          select case(ISOLVER)
          case(1) ; call solve1(         NUMBER_OF_RHS)
          case(2) ; call mumps_sc('G')
!!          case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
!!         case(3) ; call uhm_solve
!!                    call uhm_solver_flush(UHM_SOLVER_PTR)
          endselect
!
!         create tag
          itag0 = ITAG*10000 + 1*1000 + NPX*100 + NPY*10 + NPZ*1
!          
!         compute error
          call compute_error(iphys,itag0)
!
        enddo
        enddo
        enddo
!        
!       H I G H E R    O R D E R    P O L Y N O M I A L S
!
!       H1 , H(curl) , H(div)
        do ispace=1,3
!        
!       force "exact" to return appropriate higher order polynomials
        NEDELEC=ispace
!        
!       set space-specific parameters
        select case(ispace)
        case(1) ; p=0 ; ncomp=1
        case(2) ; p=1 ; ncomp=3
        case(3) ; p=1 ; ncomp=1
        endselect
!        
        do NPX=0,IORDER_TETR-p
        do NPY=0,IORDER_TETR-p
        do NPZ=0,IORDER_TETR-p
!        
        if (NPX+NPY+NPZ /= IORDER_TETR-p)  cycle
!        
          do ICOMP=1,ncomp
!  
!           solve
            select case(ISOLVER)
            case(1) ; call solve1(         NUMBER_OF_RHS)
            case(2) ; call mumps_sc('G')
!!            case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
!!            case(3) ; call uhm_solve
!!                      call uhm_solver_flush(UHM_SOLVER_PTR)
            endselect
!
!           create tag
            itag0 = ITAG*10000 + 2*1000 + NPX*100 + NPY*10 + NPZ*1
!            
!           compute error
            call compute_error(iphys,itag0)
!
          enddo
        enddo
        enddo
        enddo
!        
        enddo
!
!
!----------------------------------------------------------------------------------
!     HEXAHEDRON                                                                  |
!----------------------------------------------------------------------------------
      case('hexa')
!
!       L O W E R    O R D E R    P O L Y N O M I A L S
!
!       force "exact" to return lower order polynomials
        NEDELEC=0 ; ICOMP=0
!
        do NPX=0,IORDER_BRIC-1
        do NPY=0,IORDER_BRIC-1
        do NPZ=0,IORDER_BRIC-1
!          
!         solve
          select case(ISOLVER)
          case(1) ; call solve1(         NUMBER_OF_RHS)
          case(2) ; call mumps_sc('G')
!!          case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
!!          case(3) ; call uhm_solve
!!                    call uhm_solver_flush(UHM_SOLVER_PTR)
          endselect
!
!         create tag
          itag0 = ITAG*10000 + 1*1000 + NPX*100 + NPY*10 + NPZ*1
!          
!         compute error
          call compute_error(iphys,itag0)
!
        enddo
        enddo
        enddo
!
!       H I G H E R    O R D E R    P O L Y N O M I A L S
!
!       H1 , H(curl) , H(div)
        do ispace=1,3
!        
!       force "exact" to return appropriate higher order polynomials
        NEDELEC=ispace
!        
!       set space-specific parameters
        select case(ispace)
        case(1) ; p=0 ; ncomp=1
        case(2) ; p=1 ; ncomp=3
        case(3) ; p=1 ; ncomp=1
        endselect
!
        do NPX=0,IORDER_BRIC-p
        do NPY=0,IORDER_BRIC-p
        do NPZ=0,IORDER_BRIC-p
!          
        if ( (NPX /= IORDER_BRIC-p) .and. &
             (NPY /= IORDER_BRIC-p) .and. &
             (NPZ /= IORDER_BRIC-p)       )  cycle
!            
          do ICOMP=1,ncomp
!          
!           solve
            select case(ISOLVER)
            case(1) ; call solve1(         NUMBER_OF_RHS)
            case(2) ; call mumps_sc('G')
!!            case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
!!            case(3) ; call uhm_solve
!!                      call uhm_solver_flush(UHM_SOLVER_PTR)
            endselect
!
!           create tag
            itag0 = ITAG*10000 + 2*1000 + NPX*100 + NPY*10 + NPZ*1
!            
!           compute error
            call compute_error(iphys,itag0)
!
          enddo
        enddo
        enddo
        enddo
!
        enddo
!        
!
!----------------------------------------------------------------------------------
!     PRISM                                                                       |
!----------------------------------------------------------------------------------
      case('pris')
!
!       L O W E R    O R D E R    P O L Y N O M I A L S
!
!       force "exact" to return lower order polynomials
        NEDELEC=0 ; ICOMP=0
!
        do NPX=0,IORDER_PRIS-1
        do NPY=0,IORDER_PRIS-1
        do NPZ=0,IORDER_PRIS-1
!        
          if (NPX+NPY > IORDER_PRIS-1)  cycle
!  
!         solve
          select case(ISOLVER)
          case(1) ; call solve1(         NUMBER_OF_RHS)
          case(2) ; call mumps_sc('G')
!!          case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
!!          case(3) ; call uhm_solve
!!                    call uhm_solver_flush(UHM_SOLVER_PTR)
          endselect
!
!         create tag
          itag0 = ITAG*10000 + 1*1000 + NPX*100 + NPY*10 + NPZ*1
!          
!         compute error
          call compute_error(iphys,itag0)
!
        enddo
        enddo
        enddo
!        
!       H I G H E R    O R D E R    P O L Y N O M I A L S
!
!       H1 , H(curl) , H(div)
        do ispace=1,3
!        
!       force "exact" to return appropriate higher order polynomials
        NEDELEC=ispace
!        
!       space-specific parameters
        select case(ispace)
        case(1) ; p=0 ; ncomp=1
        case(2) ; p=1 ; ncomp=3
        case(3) ; p=1 ; ncomp=1
        endselect
!
        do NPX=0,IORDER_PRIS-p
        do NPY=0,IORDER_PRIS-p
        do NPZ=0,IORDER_PRIS-p
!        
          if ( NPX+NPY > IORDER_PRIS-p )  cycle

          if ( (NPX+NPY /= IORDER_PRIS-p) .and. &
               (NPZ     /= IORDER_PRIS-p)       )  cycle
!
          do ICOMP=1,ncomp
!  
!           solve
            select case(ISOLVER)
            case(1) ; call solve1(         NUMBER_OF_RHS)
            case(2) ; call mumps_sc('G')
!!            case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
!!            case(3) ; call uhm_solve
!!                      call uhm_solver_flush(UHM_SOLVER_PTR)
            endselect
!
!           create tag
            itag0 = ITAG*10000 + 2*1000 + NPX*100 + NPY*10 + NPZ*1
!            
!           compute error
            call compute_error(iphys,itag0)
!
          enddo
        enddo
        enddo
        enddo
!        
        enddo
!
!
!----------------------------------------------------------------------------------
!     PYRAMID                                                                     |
!----------------------------------------------------------------------------------
      case('pyra')
!
! P. Gatto , Jan 15 : this assumes that pyramid space coincides with tet space
!
!       L O W E R    O R D E R    P O L Y N O M I A L S
!
!       force "exact" to return lower order polynomials
        NEDELEC=0 ; ICOMP=0
!        
        do NPX=0,IORDER_PYRA-1
        do NPY=0,IORDER_PYRA-1
        do NPZ=0,IORDER_PYRA-1
!        
          if (NPX+NPY+NPZ > IORDER_PYRA-1)  cycle
!  
!         solve
          select case(ISOLVER)
          case(1) ; call solve1(         NUMBER_OF_RHS)
          case(2) ; call mumps_sc('G')
!!          case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
!!          case(3) ; call uhm_solve
!!                    call uhm_solver_flush(UHM_SOLVER_PTR)
          endselect
!
!         create tag
          itag0 = ITAG*10000 + 1*1000 + NPX*100 + NPY*10 + NPZ*1
!          
!         compute error
          call compute_error(iphys,itag0)
!
        enddo
        enddo
        enddo
!        
!       H I G H E R    O R D E R    P O L Y N O M I A L S
!
!       H1 , H(curl) , H(div)
        do ispace=1,3
!        
!       force "exact" to return appropriate higher order polynomials
        NEDELEC=ispace
!        
!       set space-specific parameters
        select case(ispace)
        case(1) ; p=0 ; ncomp=1
        case(2) ; p=1 ; ncomp=3
        case(3) ; p=1 ; ncomp=1
        endselect
!        
        do NPX=0,IORDER_PYRA-p
        do NPY=0,IORDER_PYRA-p
        do NPZ=0,IORDER_PYRA-p
!        
        if (NPX+NPY+NPZ /= IORDER_PYRA-p)  cycle
!        
          do ICOMP=1,ncomp
!  
!           solve
            select case(ISOLVER)
            case(1) ; call solve1(         NUMBER_OF_RHS)
            case(2) ; call mumps_sc('G')
!!            case(2) ; call mumps_solve_seq(NUMBER_OF_RHS)
!!            case(3) ; call uhm_solve
!!                      call uhm_solver_flush(UHM_SOLVER_PTR)
            endselect
!
!           create tag
            itag0 = ITAG*10000 + 2*1000 + NPX*100 + NPY*10 + NPZ*1
!            
!           compute error
            call compute_error(iphys,itag0)
!
          enddo
        enddo
        enddo
        enddo
!        
        enddo
!
      case default
        write(*,*)'test_1_aux: unknown problem!'
        stop
      endselect
!            
!
endsubroutine test_1_aux
