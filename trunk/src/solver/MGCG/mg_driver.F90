!-----------------------------------------------------------------------
!
!    routine name       - mgcg_driver
!
!-----------------------------------------------------------------------
!
!    latest revision    - Jan 2018
!
!    purpose            - debugging driver for prolongation module
!
!    arguments
!           in
!              Ref_type   - 1: uniform, 2: adaptive
!              Ref_flag   - 1: adaptive h ref, 2: adaptive p ref, 3: adaptive hp ref
!              Ref_fact   - refinement factor
!              Mg_csolve  - coarse solver selection:
!                         - 0-NO coarse solve, 1-MUMPS, 2-PARDISO
!              Mg_nrgrids - number of max grids for the MG
!              Mg_nriter  - Max number of CG itererations
!              Mg_nrsmth  - number of pre and post smoothing steps
!              Mg_tol     - MGCG stopping tolerance
!              Mg_theta   - smoother relaxation parameter
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   subroutine mg_driver(Ref_type,Ref_flag,Ref_fact,         &
                        Mg_csolve,Mg_nrgrids,Mg_nriter,Mg_nrsmth,Mg_tol,Mg_theta, &
                        Mg_export, Mg_iparAttr)
!
   use mg_data_structure, only: mg_init, mg_finalize, NRGRIDS, MAXGEN_PR, &
                                COARSE_SOLVER, MUMPS_SOLVER, PARDISO_SOLVER
   use physics,           only: NR_PHYSA
   use data_structure3D,  only: NRELES
   use constraints_info,  only: allocate_constr, store_constraints, deallocate_constr
   use patch_info,        only: patch_nodes
   use mem_alloc,         only: mg_dealloc
   use refinements,       only: enable_iso_only

!
   implicit none

   integer, intent(in)    :: Ref_type,Ref_flag,Mg_csolve,Mg_nrgrids,Mg_nrsmth, Mg_nriter
   integer, intent(in)    :: Mg_iparAttr(NR_PHYSA)
   real*8,  intent(in)    :: Ref_fact,Mg_tol,Mg_theta
   logical, intent(in)    :: Mg_export
   integer                :: igrid,nstop, iel, mdle, i, idec_sch
   integer, save          :: imat=1


!
!-----------------------------------------------------------------------
!
!..mg solver works only with isotropic refinements
!..so just to make sure...
   call enable_iso_only
!
!..set mg_parameters
   call mg_set_param(Mg_csolve,Mg_nrgrids,Mg_nrsmth,Mg_nriter,Mg_tol,Mg_theta)

!..initialize multigrid solver
   call mg_init
!
!..save mdle nodes of the coarse grid
   igrid = 1
   call fill_mdle_list(igrid)
!
!..mark current mesh
   call mark_masters(igrid)
!
   write(*,*) 'Solve coarse grid system...'

   if (COARSE_SOLVER .eq. MUMPS_SOLVER) then
      call coarse_solve_mumps('H')
   else
      call coarse_solve_pardiso
   endif
   call compute_sol_dof(igrid)
!
   do igrid = 1,NRGRIDS
      MAXGEN_PR = 1
!
      if (Mg_export) then
         write(*,*) 'Dumping mesh to matlab'
         call matlab_dump_mesh(imat)
         imat=imat+1
         write(*,*) 'Dumping solution to paraview'
         call my_paraview_driver(Mg_iparAttr)
      endif
!
      write(*,*) 'Calling refine_DPG...'
      call refine_DPG(Ref_type,Ref_flag,Ref_fact, nstop)
!
      if (igrid .lt. nrgrids) then
         write(*,*) 'Calling fill_mdle_list...'
         call fill_mdle_list(igrid+1)
      endif
!
      write(*,*) 'Create patch data structure...'
      call patch_nodes(igrid)
!
      call allocate_constr(igrid)
      write(*,*) 'Store prolongation coefficients...'
      call store_constraints(igrid)
!
      if (igrid .lt. NRGRIDS) then
         idec_sch = 1
      else
         idec_sch = 0
      endif

      write(*,*) 'Macro assembly...'
      call macro_assembly(igrid,idec_sch)
!
      write(*,*) 'patch_assembly'
      call patch_assembly(Igrid)

      write(*,*) 'pcg_solve'
      call pcg_solve(Igrid)



      if (igrid .lt. nrgrids) then
         write(*,*) 'Computing vertex patch'
         call compute_vert_patch(igrid)
!
!     ...mark this mesh as an intermediate grid
         write(*,*) 'Mark current mesh as a master'
         call mark_masters(igrid+1)
      endif
!
   enddo
!
   call mg_dealloc
!
!..finalize multigrid solver
   write(*,*) 'MG finalize'
   call mg_finalize
!
   if (Mg_export) then
      write(*,*) 'Dumping mesh to matlab'
      call matlab_dump_mesh(imat)
      imat=imat+1
      write(*,*) 'Dumping solution to paraview'
      call my_paraview_driver(Mg_iparAttr)
   endif

!
   write(*,*) 'Calling refine_DPG'
   call refine_DPG(Ref_type,Ref_flag,Ref_fact, nstop)
!
   end subroutine mg_driver
!
!
!
   subroutine mg_set_param(Mg_csolve,Mg_nrgrids,Mg_nrsmth,Mg_nriter,Mg_tol,Mg_theta)
!
   use mg_data_structure, only: NRGRIDS,   &
                                COARSE_SOLVER, MUMPS_SOLVER, PARDISO_SOLVER,NO_CSOLVE
   use pcg_info,          only: TOL, THETA, MAX_ITER, MAX_SM_ITER
!
   implicit none
!
   integer, intent(in) :: Mg_csolve, Mg_nrgrids,Mg_nrsmth,Mg_nriter
   real*8,  intent(in) :: Mg_tol, Mg_theta
!
   COARSE_SOLVER = Mg_csolve
   select case(Mg_csolve)
   case(NO_CSOLVE)
      write(*,1000) 'No coarse grid solve'
   case(MUMPS_SOLVER)
      write(*,1000) 'Coarse solver choice: MUMPS'
   case(PARDISO_SOLVER)
      write(*,1000) 'Coarse solver choice: PARDISO'
   case default
      write(*,1000) 'Wrong choice of coarse solver. Default choice: PARDISO'
      COARSE_SOLVER = PARDISO_SOLVER
   end select
!
   if (Mg_nrgrids .gt. 1) then
      write(*,1010) 'Setting max NRGRIDS     = ',Mg_nrgrids
      NRGRIDS = Mg_nrgrids
   else
      write(*,1010) 'Invalid value for Mg_nrgrids. Default choice: ', 5
      NRGRIDS = 5
   endif
!
   if (Mg_nriter .gt. 0 .and. Mg_nriter .lt. 201) then
      write(*,1010) 'Setting MAX_ITER        = ',Mg_nriter
      MAX_ITER = Mg_nriter
   else
      write(*,1010) 'Invalid value for Mg_nriter. Default choice: ', 100
      MAX_ITER = 100
   endif
!
   if (Mg_nrsmth .gt. 0 .and. Mg_nrsmth .lt. 50) then
      write(*,1010) 'Setting MAX_SM_ITER     = ',Mg_nrsmth
      MAX_SM_ITER = Mg_nrsmth
   else
      write(*,1010) 'Invalid value for Mg_nrsmth. Default choice: ', 10
      MAX_SM_ITER = 10
   endif
!
   if (Mg_tol .gt. 1d-14 .and. Mg_tol .lt. 1d-1 ) then
      write(*,1020) 'Setting TOL             = ',Mg_tol
      TOL = Mg_tol
   else
      write(*,1020) 'Invalid value for Mg_tol. Default choice: ', 1d-3
      TOL = 1d-3
   endif
!
   if (Mg_theta .gt. 0.01d0 .and. Mg_tol .lt. 0.99d0) then
      write(*,1020) 'Setting THETA           = ',Mg_theta
      THETA = Mg_theta
   else
      write(*,1020) 'Invalid value for Mg_theta. Default choice: ', 0.2d0
      THETA = 0.2d0
   endif
!
 1000 format(' mg_set_param: ',A)
 1010 format(' mg_set_param: ',A,I4)
 1020 format(' mg_set_param: ',A,ES10.3)
!
   end subroutine mg_set_param
