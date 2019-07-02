!----------------------------------------------------------------------
!                                                                     
!     program name      - main
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - main driver for the UW acoustics
!                                                                    
!----------------------------------------------------------------------
!    
   program main
!
   use environment
   use common_prob_data
   use parameters, only : NSTD_OUT, MAXbrickH
   use data_structure3D , only : NRELES, NRDOFSH, NODES
   use uhm2
   use m_assembly, ONLY: NRDOF_TOT, NRDOF_CON
!
   implicit none
!   
   real*8  :: t,greedy,rvoid,factor
   integer :: i,ic,idec,iter,nvoid,istep,nsteps,nreflag,nstop,idec_solve,j,k
   logical :: solved
   integer :: iParAttr(1), iref



!-----------------------------------------------------------------------
!                             INITIALIZATION
!-----------------------------------------------------------------------
!
!..Set common hp3d environment parameters (reads in options arguments)
   call begin_environment  ! <-- found inside src/modules/environment.F90
!
!..Set environment parameters specific to LINEAR_ELASTICITY
!..This sets default values for: FILE_REFINE,FILE_VIS,VLEVEL,
!                                FILE_CONTROL,FILE_GEOM,FILE_ERR,
!                                FILE_HISTORY,FILE_PHYS
   call set_environment  ! <-- found inside ../common/set_environment.F90
!
!..Exit if this is a "dry run".
   call end_environment  ! <-- found inside src/modules/environment.F90
!
!..print fancy header
   write(*,*)'                      '
   write(*,*)'// --  STANDARD FEM ACOUSTICS  -- //'
   write(*,*)'                      '
!
!..Initialize common library (set common parameters, load solvers, 
!                                          and create initial mesh)
   call initialize  ! <-- found inside ../common/initialize.F90
!   
!-----------------------------------------------------------------------
!                            INTERACTIVE MODE
!-----------------------------------------------------------------------
!
   write(*,*) 'OMEGA = ', OMEGA

!..display menu
   solved = .FALSE.
   idec=1 ; iref = 0
!..display menu in infinite loop
   do while(idec /= 0)
!
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      write(*,*) 'SELECT'
      write(*,*) 'QUIT ...................................0'
      write(*,*) '                                         '
      write(*,*) 'Geometry graphics (X11) ................1'
      write(*,*) 'HP3D graphics (X11) ....................2'
      write(*,*) 'Paraview ...............................3'
      write(*,*) '                                         '
      write(*,*) 'Print Data Structure arrays ...........10'
      write(*,*) 'Dumpout Data Structure ................11'
      write(*,*) '                                         '
      write(*,*) '         ---- Refinements ----           '
      write(*,*) 'Single Uniform P-refinement............20'
      write(*,*) 'Single Uniform H-refinement ...........21'
      write(*,*) 'Multi-step Uniform H-refinement .......22'
      write(*,*) '                                         '
      write(*,*) '          ----  Solvers  -----           '
      write(*,*) 'FRONTAL SOLVE PROBLEM .................30'
      write(*,*) 'MUMPS SOLVE (SEQ.) ....................40'
      write(*,*) 'MUMPS SOLVE (OMP) .....................45'
      write(*,*) 'PARDISO SOLVE (OMP)....................50'
      write(*,*) '                                         '
      write(*,*) 'EXACT ERROR ...........................60'
      write(*,*) 'MATLAB MESH PLOT.......................72'

      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      read( *,*) idec
!
      select case(idec)
!  ...QUIT
      case( 0)
         call finalize ! <-- found inside ../common/finalize.F90
      stop
!
!-----------------------------------------------------------------------
!                                 GRAPHICS
!-----------------------------------------------------------------------
!
!  ...GMP x11 graphics
      case(1) ; call graphg
!
!  ...hp3d x11 graphics
      case(2) ; call graphb
!
!  ...Paraview graphics
      case(3) 
         iParAttr = 1
         call paraview_driver(iParAttr)
!         
!-----------------------------------------------------------------------
!                                 DATA STRUCTURE
!-----------------------------------------------------------------------
!
!  ...print data structure
      case(10) ; call result
!
!  ...dump out
      case(11) ; call dumpout
!         
!-----------------------------------------------------------------------
!                               REFINEMENTS
!-----------------------------------------------------------------------
!
      case(20)
         if (.not.solved) write(*,*) 'You have not solved.'
         if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
         call global_pref
         call close
!     ...enforce max rule         
         call enforce_max_rule
!
!     ...update geometry and Dirichlet flux dof after the refinements
         call update_gdof
         call update_Ddof
!         
         solved=.FALSE.
         

!  ...Single uniform refinement
      case(21)
         if (.not.solved) write(*,*) 'You have not solved.'
         if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
         call uniform_href(IUNIFORM,1,0.25d0, nstop)
         solved=.FALSE.
!         
!  ...Multi-step uniform h refinement
       case(22)
         nsteps=0
         do while (nsteps.le.0)
            write(*,*) 'NUMBER OF REFINEMENTS>0'
            write(*,*) 'SOLVER:frontal=1,MUMPS=2,PARDISO=3,OMP MUMPS=4'
            write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
            read(*,*) nsteps,idec_solve
         enddo
         do i=0,nsteps
!        ...solve first if needed
            if (.not.solved) then
               select case(idec_solve)
               case(1)
                  call solve1(NR_RHS_PROB)
               case(2)
                  call mumps_solve_seq(NR_RHS_PROB)
               case(3)
                  call pardiso_sc
               case(4)
                  call mumps_sc('G')
               end select
            endif
!        ...say it has solved and save results to paraview file
            solved=.TRUE.
!        ...display error and refine if necessary
            if (i.ne.nsteps) then
               call uniform_href(IUNIFORM,1,0.25d0, nstop)
               if (nstop.eq.1) then
                  write(*,*) 'No elements were refined.'
                  write(*,7000) i
 7000             format('Exiting loop after ',i2,' refinements...')
                  cycle
               else
                  solved=.FALSE.
               endif
            else ! Last step only display (no refinement)
               call uniform_href(INOREFINEMENT,1,0.25d0, nstop)
            endif
        enddo
!
!-----------------------------------------------------------------------
!                                 SOLVERS
!-----------------------------------------------------------------------
!  ...Frontal Solver
      case(30)
         i=1
         if (solved) then
            write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
            read(*,*) i
         endif
         if (i.eq.1) then
            call solve1(NR_RHS_PROB)
!        ...say it has solved and save results to paraview file
            solved=.TRUE.
            ! call paraview_driver
         endif
!
!  ...MUMPS
      case(40)
         i=1
         if (solved) then
            write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
            read(*,*) i
         endif
         if (i.eq.1) then
            call uhm_time_in
            call mumps_solve_seq(NR_RHS_PROB)
            call uhm_time_out(t)
            write(*,*) 'time mumps = ', t
!        ...say it has solved and save results to paraview file
            solved=.TRUE.
        ! call paraview_driver
         endif
!
!  ...OMP MUMPS
      case(45)
         i=1
         if (solved) then
            write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
            read(*,*) i
         endif
         if (i.eq.1) then
            call mumps_sc('G')
!        ...say it has solved and save results to paraview file
            solved=.TRUE.
        ! call paraview_driver
         endif
!         
!     ...PARDISO
      case(50)
         i=1
         if (solved) then
            write(*,*) 'System already solved. Proceed anyway (Yes=1,No=0)?'
            read(*,*) i
         endif
         if (i.eq.1) then
            call pardiso_sc
!        ...say it has solved and save results to paraview file
            solved=.TRUE.
        ! call paraview_driver
         endif
!
!-----------------------------------------------------------------------
!                                 ERROR
!-----------------------------------------------------------------------
!
!    ...Compute errors of different variables
      case(60)
         call exact_error
!
      case(72)
!
         call matlab_dump_mesh(0)

      end select
!
!  ...end infinite loop
   enddo
!
!..finalize library
   call finalize ! <-- found inside ../common/finalize.F90
!
!
   end program main


!..dump mesh to a txt file for MATLAB

   subroutine matlab_dump_mesh(k)
!
   use parameters, only       : MAXbrickH
   use data_structure3D, only : NRELES, NODES

   implicit none
   integer, intent(in) :: k
   integer :: mdle, nordh,nordv,nord,i
   real*8, allocatable  :: xnod(:,:)
   character(len=20) :: str

!---------------------------------------------------------------------------------

   allocate(xnod(3,MAXbrickH)) ; xnod = 0.0d0
!   
   open(12, file='output/MATLAB/mesh_'//trim(str(k))//'.txt',       &
   status='replace',form='formatted',access='sequential')
   mdle=0;
   do i=1,NRELES
      call nelcon(mdle, mdle);
      call nodcor_vert(Mdle, xnod)
      write(12,*)                                                     &              
      (/xnod(1,1),xnod(2,1),xnod(3,1), xnod(1,2),xnod(2,2),xnod(3,2), &
        xnod(1,3),xnod(2,3),xnod(3,3), xnod(1,4),xnod(2,4),xnod(3,4), &     
        xnod(1,5),xnod(2,5),xnod(3,5), xnod(1,6),xnod(2,6),xnod(3,6), &
        xnod(1,7),xnod(2,7),xnod(3,7), xnod(1,8),xnod(2,8),xnod(3,8)/)
   enddo
!   
!..dump out the order of approximation
   mdle = 0 
   do i = 1, NRELES
      call nelcon(mdle,mdle)
      nord = NODES(mdle)%order
      call decode(nord, nordh,nordv)
      write(12,*) nordv 
   enddo
   write(12,*) NRELES
   close(12)   
   deallocate(xnod)
!
!
   end subroutine matlab_dump_mesh



!..Convert an integer to string

   character(len=20) function str(k)
! 
   integer, intent(in) :: k
   write (str, *) k
   str = adjustl(str)
! 
   end function str