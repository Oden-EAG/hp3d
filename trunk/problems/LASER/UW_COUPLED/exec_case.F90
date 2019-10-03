!----------------------------------------------------------------------
! exec_case
!----------------------------------------------------------------------
subroutine exec_case(idec)
!
   use commonParam
   use data_structure3D
   use par_mesh
   use zoltan_wrapper, only: zoltan_w_partition,zoltan_w_eval
   use mpi_param     , only: RANK,ROOT
   use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
!
   implicit none
!
   integer, intent(in) :: idec
!
   integer :: nstop
   logical :: solved
!
   integer :: flag(6)
   integer :: physNick
!
   integer :: fld,numPts
!
   integer :: src,count,ierr
!
!----------------------------------------------------------------------
!
   physNick = 1
   flag=0; flag(5)=1
!
   solved = .false.
!
   select case(idec)
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
         call href_solve(flag,physNick, nstop)
!
!  ...single anisotropic h-refinement (in z)
      case(23)
         write(*,*) 'global anisotropic h-refinement...'
         call global_href_aniso(0,1)
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
      case(50)
         write(*,*) 'computing error and residual...'
         call exact_error(flag,physNick)
!
      case(60)
         if (RANK .eq. ROOT) then
            write(*,*) 'Computing power..'
            write(*,*) 'Select field:  0 = pump'
            write(*,*) '               1 = signal'
            write(*,*) '               2 = signal and pump'
            read(*,*) fld
            write(*,*) 'Choose number sample points: 0 = default'
            read(*,*) numPts
         endif
         count = 1; src = ROOT
         call MPI_BCAST (fld   ,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
         call MPI_BCAST (numPts,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
         call get_power(fld,numPts,-1)
      case(61)
         if (RANK .eq. ROOT) then
            write(*,*) 'Computing temperature..'
            write(*,*) 'Choose number sample points: 0 = default'
            read(*,*) numPts
         endif
         count = 1; src = ROOT
         call MPI_BCAST (numPts,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
         call get_avgTemp(numPts,-1)
!
      case default
         write(*,*) 'exec_case: unknown case...'
   end select
!
end subroutine exec_case
