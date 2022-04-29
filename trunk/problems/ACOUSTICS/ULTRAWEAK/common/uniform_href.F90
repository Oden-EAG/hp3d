!--------------------------------------------------------------------
!
!     routine name      - uniform_href
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!     purpose:          - performs single uniform h-refinement
!
!     arguments:
!        in:
!             Irefine   = 0 if no refinement -> only display
!                       = 1 if uniform refinements
!             Nreflag   = 1 for uniform isotropic h-refinements
!                       = 2 for uniform anisotropic h-refinements (in z)
!             Factor    - if element error \ge Factor*max_error
!                         then the element is refined
!
!---------------------------------------------------------------------
subroutine uniform_href(Irefine,Nreflag,Factor)
!
   use control
   use data_structure3D
   use environment
   use common_prob_data
   use assembly_sc, only: NRDOF_CON,NRDOF_TOT
   use par_mesh   , only: DISTRIBUTED,HOST_MESH
   use mpi_param  , only: ROOT,RANK
   use MPI        , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD,MPI_REAL8
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real(8),  intent(in) :: Factor
!
   integer, parameter :: max_step = 20
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real(8), dimension(max_step), save :: rate_mesh
   real(8), dimension(max_step), save :: error_mesh
   real(8), dimension(max_step), save :: rel_error_mesh
   real(8), dimension(max_step), save :: rate_error_mesh
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
!
   integer :: iflag(1)
   real(8) :: errorH,errorE,errorV,errorQ
   real(8) :: rnormH,rnormE,rnormV,rnormQ
   real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd
   integer :: i,iel,mdle,count,ierr,nrelem_ref,kref
   integer :: iprint
!
   real(8) :: MPI_Wtime,start_time,end_time
!
!-----------------------------------------------------------------------
!
!..initialize
   iprint = 0
!
!..increase step if necessary
!..irefineold=0 means no refinement was made in the previous step
!..if first call or if a refinement was made, increase step
   if ((istep.eq.0).or.(irefineold.ne.0)) istep=istep+1
   irefineold=Irefine
!
!..get dof count from solver
   nrdof_tot_mesh(istep) = NRDOF_TOT
   nrdof_con_mesh(istep) = NRDOF_CON
   if (NRDOF_TOT .eq. 0) then
      istep=istep-1
      goto 90
   endif
!
!..field variables flag
   iflag(1) = 1
!
!..fetch active elements
   if (.not. DISTRIBUTED) then
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
   error_subd = 0.d0; rnorm_subd = 0.d0
!
!$OMP PARALLEL DEFAULT(PRIVATE)              &
!$OMP SHARED(NRELES_SUBD,ELEM_SUBD,iflag)    &
!$OMP REDUCTION(+:error_subd,rnorm_subd)
!$OMP DO
   do iel=1,NRELES_SUBD
      call element_error(ELEM_SUBD(iel),iflag,           &
                         errorH,errorE,errorV,errorQ,    &
                         rnormH,rnormE,rnormV,rnormQ)
      error_subd = error_subd + errorQ
      rnorm_subd = rnorm_subd + rnormQ
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
   error_tot = 0.d0; rnorm_tot = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_REDUCE(error_subd,error_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rnorm_subd,rnorm_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
   else
      error_tot = error_subd
      rnorm_tot = rnorm_subd
   endif
!
   if (RANK .ne. ROOT) goto 90
!
!..update and display convergence history
   if (NEXACT.ge.1) then
      error_mesh(istep) = sqrt(error_tot)
      rel_error_mesh(istep) = sqrt(error_tot/rnorm_tot)
   endif
!
!..compute decrease rate for the residual and error
   select case(istep)
   case(1)
      rate_mesh(istep) = 0.d0
      if (NEXACT.ge.1) rate_error_mesh(istep) = 0.d0
   case default
      if (NEXACT.ge.1) then
         rate_error_mesh(istep) = &
         log(rel_error_mesh(istep-1)/rel_error_mesh(istep))/  &
         log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      endif
   end select
!
!..print out the history of refinements
   if (NEXACT.gt.0) then
      write(*,*)
      write(*,*) 'HISTORY OF REFINEMENTS'
      write(*,110)
  110 format(' mesh |',' nrdof_tot |',' nrdof_con | ',' |',   &
             ' field error  |','rel field error|','   rate ')
      write(*,*)
!
      do i=1,istep
         write(*,120) i, nrdof_tot_mesh(i), nrdof_con_mesh(i), &
            error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
  120    format(2x,i2,'  | ',2(i9,' | '), 2(' | ',es12.5),'  |',f7.2)
         if (i .eq. istep) write(*,*)
      enddo
   endif
!
 90 continue
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
  2010 format(A,I3,A)
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
         if (Nreflag .eq. 1) then
!        ...isotropic href
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
            call global_href
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
            call update_gdof
            call update_Ddof
         elseif (Nreflag .eq. 2) then
!        ...anisotropic href (z)
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
            call global_href_aniso(0,1) ! refine in z
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
            call update_gdof
            call update_Ddof
         elseif (Nreflag .eq. 3) then
!        ...anisotropic href (z)
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
            nrelem_ref = NRELES
            do iel=1,nrelem_ref
               mdle = ELEM_ORDER(iel)
               select case(NODES(mdle)%type)
                  case('mdlb'); kref = 1 ! refine in z
                  case('mdlp'); kref = 1 ! refine in z
                  case default
                     write(*,*) 'href WARNING: anisoref for unexpected type. stop.'
                     stop 1
               end select
               call refine(mdle,kref)
            enddo
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2030) ' # of current elements, nodes = ', NRELES, NRNODS
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
            call close_mesh
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
               write(*,2025) end_time-start_time
  2020         format(' global ref : ',f12.5,'  seconds')
  2025         format(' close mesh : ',f12.5,'  seconds')
  2030         format(A,I8,', ',I9)
            endif
            call update_gdof
            call update_Ddof
         else
            write(*,*) 'href WARNING: unexpected Nreflag'
         endif
   end select
!
end subroutine uniform_href
!
!
subroutine href_solve()
!
   use common_prob_data
   use MPI           , only: MPI_COMM_WORLD,MPI_INTEGER
   use mpi_param     , only: RANK,ROOT
   use par_mesh      , only: DISTRIBUTED,HOST_MESH,distr_mesh
   use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_partition
!
   implicit none
!
   integer :: nsteps,count,src,ierr,i
!
   if (RANK .eq. ROOT) then
      nsteps=0
      do while (nsteps.le.0)
         write(*,*) 'Provide: number of uniform h-refinements'
         read(*,*) nsteps
      enddo
   else
      write(6,405) '[', RANK, '] : ','Waiting for broadcast from master...'
 405  format(A,I4,A,A)
   endif
   count = 1; src = ROOT
   call MPI_BCAST (nsteps,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
!
!..initial solve
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      call par_mumps_sc('G')
   else
      call mumps_sc('G')
   endif
!
!..do refinements and solve
   do i=0,nsteps
!  ...display error and refine if necessary
      if (i.ne.nsteps) then
!     ...Uniform refinement and solve
         call uniform_href(IUNIFORM,1,0.25d0)
         if (DISTRIBUTED .and. (.not. HOST_MESH)) then
            !call zoltan_w_set_lb(1)
            !call distr_mesh
            !call print_partition
            call par_mumps_sc('G')
         else
            call mumps_sc('G')
         endif
      else
!     ...Last step only display (no refinement)
         call uniform_href(INOREFINEMENT,1,0.25d0)
      endif
   enddo
!
end subroutine href_solve

