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
!             Nreflag   = 1 for h-refinements
!                       = 2 for p-refinements
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
   use MPI        , only: MPI_SUM,MPI_COMM_WORLD,MPI_REAL8
!
   implicit none
!   
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real*8,  intent(in)  :: Factor
!
   integer, parameter :: max_step = 20
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real*8,  dimension(max_step), save :: rate_mesh
   real*8,  dimension(max_step), save :: error_mesh
   real*8,  dimension(max_step), save :: rel_error_mesh
   real*8,  dimension(max_step), save :: rate_error_mesh
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
!
   integer :: iflag(1)
   integer :: mdle_list(NRELES)
   real*8  :: errorH,errorE,errorV,errorQ
   real*8  :: rnormH,rnormE,rnormV,rnormQ
   real*8  :: error,rnorm,error_subd,rnorm_subd
   integer :: i,iel,mdle,subd,nreles_subd,count,ierr
   integer :: iprint
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
!
!..field variables flag
   iflag(1) = 1
!
!..fetch active elements
   mdle = 0
   if ((DISTRIBUTED) .and. (.not. HOST_MESH)) then
      nreles_subd = 0
      do iel=1,NRELES
         call nelcon(mdle, mdle)
         call get_subd(mdle, subd)
         if (subd .eq. RANK) then
            nreles_subd = nreles_subd + 1
            mdle_list(nreles_subd) = mdle
         endif
      enddo
   else
      if (RANK .ne. ROOT) goto 90
      do iel=1,NRELES
         call nelcon(mdle, mdle)
         mdle_list(iel) = mdle
      enddo
      nreles_subd = NRELES
   endif
!
   error_subd = 0.d0; rnorm_subd = 0.d0
!
!$OMP PARALLEL DEFAULT(PRIVATE)              &
!$OMP SHARED(nreles_subd,iflag,mdle_list)    &
!$OMP REDUCTION(+:error_subd,rnorm_subd)
!$OMP DO
   do iel=1,nreles_subd
      call element_error(mdle_list(iel),iflag,           &
                         errorH,errorE,errorV,errorQ,    &
                         rnormH,rnormE,rnormV,rnormQ)
      error_subd = error_subd + errorH
      rnorm_subd = rnorm_subd + rnormH
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
   error = 0.d0; rnorm = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_REDUCE(error_subd,error,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(rnorm_subd,rnorm,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
   else
      error = error_subd
      rnorm = rnorm_subd
   endif
!
!..update and display convergence history
   if (NEXACT.ge.1) then
      error_mesh(istep) = sqrt(error)
      rel_error_mesh(istep) = sqrt(error/rnorm)
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
   if ((NEXACT.gt.0) .and. (RANK.eq.ROOT)) then
      write(*,*)
      write(*,*) 'HISTORY OF REFINEMENTS'
      write(*,110)
  110 format(' mesh |',' nrdof_tot|',' nrdof_con| ',' |',   &
             ' field error  |','rel field error|','   rate ')
      write(*,*)
!
      do i=1,istep
         write(*,120) i, nrdof_tot_mesh(i), nrdof_con_mesh(i), &
            error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
  120    format(2x,i2,'  | ',2(i8,' | '), 2(' | ',es12.5),'  |',f7.2)
      enddo
   endif
!
 90 continue
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
         call global_href
         call close_mesh
         call update_gdof
         call update_ddof
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

