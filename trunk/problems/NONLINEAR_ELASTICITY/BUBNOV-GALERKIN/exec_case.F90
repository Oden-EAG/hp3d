!----------------------------------------------------------------------
! exec_case
!----------------------------------------------------------------------
subroutine exec_case(idec)
!
   use data_structure3D
   use par_mesh
   use mpi_param
   use mpi
   use common_prob_data
   use zoltan_wrapper, only: zoltan_w_partition,zoltan_w_eval
!
   implicit none
!
   integer, intent(in) :: idec
!
   logical :: solved , Lsflag_tmp,Bfgs_tmp
   integer :: mdle_subd(NRELES)
   integer :: i,mdle,kref,src,count,ierr

   integer :: N_incr,Nupdate,Maxiter
   real*8  :: Rel_tol,Abs_tol,Step_size
!
!----------------------------------------------------------------------
!
   solved = .false.
!
   select case(idec)
!
!  ...paraview graphics
      case(3)
         call my_paraview_driver(NR_COMP)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!  ...print data structure (interactive)
      case(10); call result
!
!  ...print general data structure info
      case(11)
         write(*,110) NRELIS,NRELES,NRNODS
 110     format(' NRELIS,NRELES,NRNODS            = ',3I10)
         write(*,111) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
 111     format(' NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ = ',4I10)
         write(*,112) MAXNODS,NPNODS
 112     format(' MAXNODS,NPNODS                  = ',2I10)
!
!  ...print current partition (elems)
      case(15)
         write(*,*) 'printing current partition (elems)...'
         call print_partition
!
!  ...print current subdomains (nodes)
      case(16)
         write(*,*) 'printing current subdomains (nodes)...'
         call print_subd
!
!  ...print current partition coordinates
      case(17)
         write(*,*) 'printing current partition coordinates...'
         call print_coord
!
!  ...single uniform h-refinement
      case(20)
         write(*,*) 'global h-refinement...'
         call global_href
         call update_gdof
         call update_Ddof
!
!  ...single uniform p-refinement
      case(21)
         write(*,*) 'global p-refinement...'
         call global_pref
         call update_gdof
         call update_Ddof
!
!  ...Multi-step uniform h refinement
      case(22)
         call href_solve
!
!  ...single anisotropic h-refinement (in z)
      case(23)
         write(*,*) 'global anisotropic h-refinement...'
         call global_href_aniso(0,1)
         call update_gdof
         call update_Ddof
!
      case(26)
         if (RANK.eq.ROOT) then
            write(*,*) 'Select on mdle node from the list: '
            do i=1,NRELES
               write(*,2610) ELEM_ORDER(i)
          2610 format(I6)
            enddo
            read(*,*) mdle
         endif
         if (NUM_PROCS .gt. 1) then
            count = 1; src = ROOT
            call MPI_BCAST (mdle,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
         endif
         kref = 111
         call refine(mdle,kref)
         call close_mesh
         call update_gdof
         call update_Ddof
!
!  ...distribute mesh
      case(30)
         write(*,*) 'distribute mesh...'
         call distr_mesh
!
!  ...collect dofs on ROOT processor
      case(31)
         write(*,*) 'collecting dofs on ROOT...'
         call collect_dofs
!
!  ...suggest new mesh partition (Zoltan)
      case(32)
         if (DISTRIBUTED) then
            write(*,*) 'computing new mesh partition (Zoltan)...'
            call zoltan_w_partition(mdle_subd)
         else
            write(*,*) 'distribute mesh first to use Zoltan...'
         endif
!
!  ...evaluate current partition
      case(33)
         if (DISTRIBUTED) then
            write(*,*) 'evaluating current partition...'
            call zoltan_w_eval
         else
            write(*,*) 'distribute mesh first to use Zoltan...'
         endif
!
!  ...run mesh verification routines
      case(35)
         write(*,*) 'verify distributed mesh consistency...'
         call par_verify
!
!  ...solve problem with omp_mumps (OpenMP MUMPS)
      case(40)
         write(*,*) 'calling MUMPS (MPI) solver...'
         call par_mumps_sc('G')
!
!  ...solve problem with par_mumps (MPI MUMPS)
      case(41)

         write(*,*) '::::::::::::::::::::::::::::::::'
         write(*,*) 'NON-LINEAR SOLVER BASED ON MUMPS'
         write(*,*) ''
         write(*,*) 'PLEASE ENTER THE FOLLOWING PARAMETERS'
         write(*,*) 'Number of load steps'
         read (*,*) N_incr
         write(*,*) 'Load step size'
         read (*,*) Step_size
         write(*,*) 'Absolute tolerance for residual'
         read (*,*) Abs_tol
         write(*,*) 'Relative tolerance for residual'
         read (*,*) Rel_tol
         write(*,*) 'Maximum number of (quasi)Newton-Raphson iterations per load step'
         read (*,*) Maxiter
         write(*,*) 'Update matrix every x iterations (1 means full N-R method)'
         read (*,*) Nupdate
         write(*,*) 'Use Line-Search (T/F)'
         read (*,*) Lsflag_tmp
         write(*,*) 'Use BFGS updates (T/F)'
         read (*,*) Bfgs_tmp
         
         call nl_elast_solve(N_incr,Step_size,Abs_tol,Rel_tol,Maxiter,Nupdate,Lsflag_tmp,Bfgs_tmp)



!  ...solve problem with par_mumps (MPI MUMPS)
      case(49)

         write(*,*) '::::::::::::::::::::::::::::::::'
         write(*,*) 'NON-LINEAR SOLVER BASED ON MUMPS'

         N_incr = 5
         Step_size = 2.d-3
         Abs_tol = 1.d-12
         Rel_tol = 1.d-8
         Maxiter = 100
         Nupdate = 10
         Lsflag_tmp = .true.
         Bfgs_tmp   = .true.
         call nl_elast_solve(N_incr,Step_size,Abs_tol,Rel_tol,Maxiter,Nupdate,Lsflag_tmp,Bfgs_tmp)



!
!  ...solve problem with pardiso (OpenMP)
      case(42)
         write(*,*) 'calling Pardiso (OpenMP) solver...'
         call pardiso_sc('G')
!
!  ...solve problem with Frontal solver (sequential)
      case(43)
         write(*,*) 'calling Frontal (Seq) solver...'
         call solve1(1)
!
!  ...solve problem with omp_mumps (OpenMP MUMPS)
      case(44)
         write(*,*) 'calling MUMPS (MPI) nested dissection solver...'
         call par_nested('G')
!
!  ...solve problem with PETSc solver (MPI)
      case(45)
         write(*,*) 'calling PETSc (MPI) solver...'
         call petsc_solve('G')
!
      case(50)
         write(*,*) 'computing error and residual...'
         call exact_error
!
      case(60)
         write(*,*) 'flushing dof'
         do i=1,NRNODS
            if(associated(NODES(i)%dof)) then
               if(associated(NODES(i)%dof%coord)) then
                  NODES(i)%dof%coord = 0.d0
               endif
               if(associated(NODES(i)%dof%zdofH)) then
                  NODES(i)%dof%zdofH = 0.d0
               endif
            endif
         enddo
         call update_gdof
         call update_Ddof
!
      case default
         write(*,*) 'exec_case: unknown case...'
   end select
!
end subroutine exec_case

