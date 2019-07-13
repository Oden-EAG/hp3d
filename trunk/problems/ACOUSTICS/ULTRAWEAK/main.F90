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
   use data_structure3D , only : NRELES, NRDOFSH, NODES, NRNODS
   use m_assembly, ONLY: NRDOF_TOT, NRDOF_CON, IPRINT_TIME, STORE_STC
   use mg_data_structure
!
   implicit none
!   
   real*8  :: t,greedy,rvoid,factor
   integer :: i,ic,idec,iter,nvoid,istep,nsteps,nreflag,nstop,idec_solve,j,k
   logical :: solved
   integer :: iParAttr(3), iref
   integer :: mdle, nord, iel, nod,nodp, nrule
!
!..initialization
!
!..Set common hp3D environment parameters (reads in options arguments)
   call begin_environment  ! <-- found inside src/modules/environment.F90
!
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
   write(*,*)'// --  ULTRAWEAK DPG FOR ACOUSTICS  -- //'
   write(*,*)'                      '
!
!..Initialize common library (set common parameters, load solvers, 
!                                          and create initial mesh)
   call initialize  ! <-- found inside ../common/initialize.F90
!   
!..interactive mode
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
      write(*,*) 'Single Adaptive H-Refinements .........22'
      write(*,*) 'Multi-step Uniform H-refinement .......23'
      write(*,*) 'Multi-step Adaptive hp-refinement......24'
      write(*,*) '                                         '
      write(*,*) '          ----  Solvers  -----           '
      write(*,*) 'FRONTAL SOLVE PROBLEM .................30'
      write(*,*) 'MUMPS SOLVE (SEQ.) ....................40'
      write(*,*) 'MUMPS SOLVE (OMP) .....................45'
      write(*,*) 'PARDISO SOLVE (HERM OMP)...............50'
      write(*,*) '                                         '
      write(*,*) 'EXACT ERROR ...........................60'
!     write(*,*) 'Write error to file ...................61'
      write(*,*) 'RESIDUAL ..............................70'
      write(*,*) 'MATLAB MESH PLOT.......................72'
      write(*,*) 'PROLONGATION DEBUG....................100'

      ! write(*,*) 'PRINT NOD COORD.......................200'
      ! write(*,*) 'FIND_FATHER...........................201'
      ! write(*,*) 'HASHING TEST...........................73'
      ! write(*,*) 'MY TESTS..............................100'

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
!  ...GRAPHICS
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
         iParAttr = (/0,0,1/)
         call paraview_driver(iParAttr)
!         
!-----------------------------------------------------------------------
!   ...DATA STRUCTURE
!-----------------------------------------------------------------------
!
!  ...print data structure
      case(10) ; call result
!
!  ...dump out
      case(11) ; call dumpout
!         
!-----------------------------------------------------------------------
!  ...REFINEMENTS
!-----------------------------------------------------------------------
!
      case(20)

        mdle=0
        do iel=1,NRELES
           call nelcon(mdle, mdle)
           nord = NODES(mdle)%order
           select case(NODES(mdle)%type)
              case('mdln','mdld'); nord = nord+1
              case('mdlp'); nord = nord+11
              case('mdlb'); nord = nord+111
           end select
           call nodmod(mdle, nord)
        enddo

!          if (.not.solved) write(*,*) 'You have not solved.'
!          if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
!          call close
! !     ...enforce max rule         
!          ! call enforce_max_rule
! !
! !     ...update geometry and Dirichlet flux dof after the refinements
!          call update_gdof
!          call update_Ddof
! !         
!          solved=.FALSE.
         

!  ...Single uniform refinement
      case(21)
         if (.not.solved) write(*,*) 'You have not solved.'
         if (.not.solved) write(*,*) 'Mesh will uniformly refine anyway.'
         call refine_DPG(IUNIFORM,1,0.25d0, nstop)
         solved=.FALSE.
!         
!     ...Single adaptive refinement
      case(22)
         if (.not.solved) then
            write(*,*) 'You have not solved. Cannot adaptively refine.'
         else
         nreflag=0
            do while ((nreflag.ne.1).and.(nreflag.ne.2).and.(nreflag.ne.3))
               write(*,*) 'REFINEMENT FLAG:h-refine=1,p-refine=2,hp-refine=3'
               write(*,*) '0.d0<FACTOR<1.d0'
               write(*,*) 'Provide: REFINEMENT FLAG, FACTOR'
               read(*,*) nreflag,factor
            enddo
            call refine_DPG(IADAPTIVE,nreflag,factor, nstop)
            if (nstop.eq.1) write(*,*) 'No elements were refined.'
            if (nstop.eq.0) solved=.FALSE.
         endif
!         
!     ...Multi-step uniform h refinement
         case(23)
            STORE_STC = .TRUE.

            nsteps=0
            do while (nsteps.le.0)
               write(*,*) 'NUMBER OF REFINEMENTS>0'
               write(*,*) 'SOLVER:frontal=1,MUMPS=2,PARDISO=3,OMP MUMPS=4'
               write(*,*) 'Provide: NUMBER OF REFINEMENTS, SOLVER'
               read(*,*) nsteps,idec_solve
            enddo
            do i=0,nsteps
!           ...solve first if needed
               if (.not.solved) then
                  select case(idec_solve)
                  case(1)
                     call solve1(NR_RHS_PROB)
                  case(2)
                     call mumps_solve_seq(NR_RHS_PROB)
                  case(3)
                     IPRINT_TIME = 1
                     ! call pardiso_sc('H')
                     call pardiso_interface
                     ! IPRINT_TIME = 0
                  case(4)
                     call mumps_sc('H')
                  end select
               endif
!           ...say it has solved and save results to paraview file
               solved=.TRUE.
!           ...display error and refine if necessary
               if (i.ne.nsteps) then
                  call refine_DPG(IUNIFORM,1,0.25d0, nstop)
                  if (nstop.eq.1) then
                     write(*,*) 'No elements were refined.'
                     write(*,7000) i
 7000                format('Exiting loop after ',i2,' refinements...')
                     cycle
                  else
                     solved=.FALSE.
                  endif
               else ! Last step only display (no refinement)
                  call refine_DPG(INOREFINEMENT,1,0.25d0, nstop)
               endif
            enddo
!        ...Multi-step adaptive hp refinement

         case(24)
            call mg_init
            nreflag=0
            do while ((nreflag.ne.1).and.(nreflag.ne.2).and.   &
                      (nreflag.ne.3).and.(nreflag.ne.4))
               write(*,*) 'NUMBER OF REFINEMENTS>0'
               write(*,*) 'REFINEMENT FLAG:'
               write(*,*) 'h-refine=1,p-refine=2,hp-refine(max)=3,hp-refine(min)=4'
               write(*,*) '0.d0<FACTOR<1.d0'
               write(*,*) 'SOLVER:frontal=1,MUMPS=2,PARDISO=3,OMP MUMPS=4 '
               write(*,*) 'Provide: NUMBER OF REFINEMENTS,',  &
                          ' REFINEMENT FLAG, FACTOR, SOLVER'
               read(*,*) nsteps,nreflag,factor,idec_solve
            enddo
            do i=0,nsteps
               iref = iref+1
!        ...solve first if needed
               if (.not.solved) then
                  select case(idec_solve)
                  case(1)
                     call solve1(NR_RHS_PROB)
                  case(2)
                     IPRINT_TIME = 1
                     call mumps_interf(NR_RHS_PROB)
                     IPRINT_TIME = 0

                  case(3)
                     IPRINT_TIME = 1
                     call pardiso_sc('H')
                     ! call pardiso_interface
                     IPRINT_TIME = 0
                  case(4)
                      IPRINT_TIME = 1
                     call mumps_sc('H')
                      IPRINT_TIME = 0
                  end select
               endif
!           ...say it has solved and save results to paraview file
               solved=.TRUE.
               j = mod(iref-1,1)
               ! if (iref-1 .gt. 1) then 
               !    if (mod(iref-1,30).eq.0) then 
               !       write(*,*) 'dumping out hp3d'
               !       write(*,*) 'iref-1, mod(iref-1,30) = ', iref-1, mod(iref-1,30)   
               !       call dumpout
               !    endif  
               ! endif   
 !               if (j .eq. 0) then 
 !                  write(*,1000) iref-1
 ! 1000             format(' Dumping to paraview: mesh = ',i3)          
 !                  iParAttr = (/0,0,1/)
 !                  call paraview_driver(iParAttr)
 !                  write(*,1001) iref-1
 ! 1001             format(' Dumping to MATLAB: mesh =  ',i3)          
 !                  call matlab_dump_mesh(iref-1) 
 !               endif  

!           ...then display error and refine if necessary
               if (i.ne.nsteps) then
                  do nod = 1, NRNODS
                     NODES_MG(nod)%orderC = NODES(nod)%order
                  enddo   
                  call refine_DPG(IADAPTIVE,nreflag,factor, nstop)
                  do nod = 1, NRNODS
                     if (NODES(nod)%order-NODES_MG(nod)%orderC .lt. 0 ) then
                        write(*,*) 'nod, type,act, orderC, order = ' , &
                        nod, NODES(nod)%type, NODES(nod)%act,    &
                              NODES_MG(nod)%orderC, NODES(nod)%order
                        call pause
                     endif
                  enddo      
                  if (nstop.eq.1) then
                     write(*,*) 'No elements were refined.'
                     write(*,7000) i
                     cycle
                  else
                     solved=.FALSE.
                  endif
               else ! Last step only display (no refinement)
                  call refine_DPG(INOREFINEMENT,nreflag,factor, nstop)
               endif
            enddo

            call mg_finalize

!
!
!-----------------------------------------------------------------------
!  ...SOLVERS
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
            !call uhm_time_in
            call mumps_solve_seq(NR_RHS_PROB)
            !call uhm_time_out(t)
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
            STORE_STC = .TRUE.
            ! call mumps_sc('H')
            call mumps_sc('H')
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
            STORE_STC = .TRUE.
            ! call pardiso_sc('H')
            call pardiso_interface
!        ...say it has solved and save results to paraview file
            solved=.TRUE.
        ! call paraview_driver
         endif
!
!-----------------------------------------------------------------------
!                                 ERROR
!-----------------------------------------------------------------------
!
!    ...Compute error
      case(60)
         call exact_error
!
!     ...Compute residual
      case(70)
         call compute_residual
!
      case(72)
!
         call matlab_dump_mesh(0)
     

      case(100)

 10      continue
         write(*,*) 'Choose max/min rule 3/4 and number of refinements'
         read(*,*) nrule, nsteps
         if (nrule .ne. 1 .and. nrule .ne. 3 .and. nrule .ne. 4) go to 10
         do i = 0,nsteps
            call debug_driver(nrule)
            ! call graphb
         enddo

      case(200)

         call find_coord

      case(201)

         do nod = 1, NRNODS
            if (NODES(nod)%act .eq. 1 .and. NODES(nod)%type .eq. 'medg') then
               call locate_father(nod,nodp)
               write(*,*) 'nod, nodp = ', nod, nodp
            endif   
         enddo   

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


!..dump mesh to a .txt file for MATLAB

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

!
!..Convert an integer to string

   character(len=20) function str(k)
! 
   integer, intent(in) :: k
   write (str, *) k
   str = adjustl(str)
! 
   end function str



   subroutine find_coord

   use data_structure3D
   implicit none
   character(len=4)  :: type
   integer           :: norientl(27),nodesl(27) 
   integer           :: mdle, nrv, nre,nrf,iel,i
!
!------------------------------------------------------------------------------
!
   mdle = 0
   do iel = 1,NRELES
      call nelcon(mdle,mdle)
      type = NODES(mdle)%type
      nrv = NVERT(Type)
      nre = NEDGE(Type)
      nrf = NFACE(Type)
!
      call elem_nodes(mdle,nodesl,norientl)
      write(*,*) 'element number, mdle = ', iel, mdle

      write(*,9998)
 9998 format ('vertex number,  vertex coord')
 
 !      do i = 1, nrv
 !         write(*,9999)  nodesl(i), floor(NODES(nodesl(i))%coord)
 ! 9999 format(i8,'--->', '(',i1,',',i1,',',i1,')')        
 !      enddo   




      ! write(*,*) 'edges    = ', nodesl(nrv+1:nrv+nre)
      ! write(*,*) 'faces    = ', nodesl(nrv+nre+1:nrv+nre+nrf)
      ! write(*,*) 'mdle     = ', nodesl(nrv+nre+nrf+1)
      call pause
   enddo



   end subroutine find_coord
