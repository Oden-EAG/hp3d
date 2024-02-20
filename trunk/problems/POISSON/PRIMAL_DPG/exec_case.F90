!----------------------------------------------------------------------
! exec_case
!----------------------------------------------------------------------
subroutine exec_case(idec)
!
   use data_structure3D
   use par_mesh
   use mpi_wrapper
   use common_prob_data
   use paraview      , only: paraview_select_attr
   use zoltan_wrapper, only: zoltan_w_partition,zoltan_w_eval
!
   implicit none
!
   integer, intent(in) :: idec
!
   logical :: solved
   integer :: mdle_subd(NRELES)
   integer :: i,mdle,kref,src,count,ierr,nord
   logical :: iPvAttr(2)
   real(8) :: res
!
!----------------------------------------------------------------------
!
   solved = .false.
!
   select case(idec)
!
!  ...paraview graphics
      case(3)
         iPvAttr(1:2) = (/.true.,.false./) ! write field output only
         call paraview_select_attr(iPvAttr)
         call paraview_driver
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
         call par_mumps_sc('H')
!
!  ...solve problem with par_mumps (MPI MUMPS)
      case(41)
         write(*,*) 'calling MUMPS (OpenMP) solver...'
         call mumps_sc('H')
!
!  ...solve problem with pardiso (OpenMP)
      case(42)
         write(*,*) 'calling Pardiso (OpenMP) solver...'
         call pardiso_sc('H')
!
!  ...solve problem with Frontal solver (sequential)
      case(43)
         write(*,*) 'calling Frontal (Seq) solver...'
         call solve1(1)
!
!  ...solve problem with omp_mumps (OpenMP MUMPS)
      case(44)
         write(*,*) 'calling MUMPS (MPI) nested dissection solver...'
         call par_nested('H')
!
!  ...solve problem with PETSc solver (MPI)
      case(45)
         write(*,*) 'calling PETSc (MPI) solver...'
         call petsc_solve('P')
!
      case(50)
         write(*,*) 'computing error...'
         call exact_error
!
      case(51)
         write(*,*) 'computing residual...'
         call residual(res)
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
      case(65)
         if (RANK.eq.ROOT) then
            write(*,*) 'Select a mdle node from the list: '
            do i=1,NRELES
               write(*,2610) ELEM_ORDER(i)
            enddo
            read(*,*) mdle
         endif
         if (NUM_PROCS .gt. 1) then
            count = 1; src = ROOT
            call MPI_BCAST (mdle,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
         endif
         if (RANK.eq.ROOT) then
            write(*,*) 'Select corresponding mdle node refinement:'
            read(*,*) nord
         endif
         if (NUM_PROCS .gt. 1) then
            count = 1; src = ROOT
            call MPI_BCAST (nord,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
         endif
         call nodmod(mdle,nord)
!
      case default
         write(*,*) 'exec_case: unknown case...'
   end select
!
end subroutine exec_case

