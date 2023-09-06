!--------------------------------------------------------------------
!
!     routine name      - uniform_href
!
!--------------------------------------------------------------------
!
!     latest revision:  - May 2023
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
   use MPI        , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD, &
                          MPI_REAL8,MPI_Wtime
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real(8), intent(in)  :: Factor
!
   integer, parameter :: max_step = 20
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real(8), dimension(max_step), save :: rate_mesh
   real(8), dimension(max_step), save :: error_mesh
   real(8), dimension(max_step), save :: rel_error_mesh
   real(8), dimension(max_step), save :: rate_error_mesh
!
   real(8), dimension(max_step), save :: residual_mesh
   real(8), dimension(max_step), save :: rate_residual_mesh
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
!
   integer :: iflag(NR_PHYSA)
   real(8) :: errorH,errorE,errorV,errorQ
   real(8) :: rnormH,rnormE,rnormV,rnormQ
   real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd
   integer :: i,iel,mdle,count,ierr,nrelem_ref,kref
   integer :: iprint
!
   real(8) :: res
   real(8) :: resid_subd,resid_tot
   real(8) :: elem_resid
   integer :: elem_ref_flag
!
   real(8) :: start_time,end_time
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
!..compute the error only for L2 field
   iflag(1:NR_PHYSA) = (/0,0,1,1/)
!
!..fetch active elements
   if (.not. DISTRIBUTED) then
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
   error_subd = 0.d0; rnorm_subd = 0.d0
   resid_subd = 0.d0
!
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(errorH,errorE,errorV,errorQ,            &
!$OMP         rnormH,rnormE,rnormV,rnormQ)            &
!$OMP SCHEDULE(DYNAMIC)                               &
!$OMP REDUCTION(+:resid_subd,error_subd,rnorm_subd) 
   do iel=1,NRELES_SUBD
      if(NEXACT .gt. 0) then 
         call element_error(ELEM_SUBD(iel),iflag,           &
                            errorH,errorE,errorV,errorQ,    &
                            rnormH,rnormE,rnormV,rnormQ)
         error_subd = error_subd + errorQ
         rnorm_subd = rnorm_subd + rnormQ
      endif

      call elem_residual(ELEM_SUBD(iel), elem_resid,elem_ref_flag)
      resid_subd = resid_subd + elem_resid
   enddo
!$OMP END PARALLEL DO
!
   error_tot = 0.d0; rnorm_tot = 0.d0;resid_tot = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      if(NEXACT .gt. 0) then 
         call MPI_REDUCE(error_subd,error_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(rnorm_subd,rnorm_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      endif

      call MPI_REDUCE(resid_subd,resid_tot,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
!
   else
      if(NEXACT .gt. 0) then 
         error_tot = error_subd
         rnorm_tot = rnorm_subd
      endif

      resid_tot = resid_subd
   endif
!
   if (RANK .ne. ROOT) goto 90
!
!..update and display convergence history
   if (NEXACT.gt.0) then
      error_mesh(istep) = sqrt(error_tot)
      rel_error_mesh(istep) = sqrt(error_tot/rnorm_tot)
   endif
   residual_mesh(istep) = sqrt(resid_tot)
!
!..compute decrease rate for the residual and error
   select case(istep)
   case(1)
!      
      rate_mesh(istep) = 0.d0
      if (NEXACT.gt.0) rate_error_mesh(istep) = 0.d0
!
      rate_residual_mesh(istep) = 0.d0
      if (NEXACT.gt.0) rate_residual_mesh(istep) = 0.d0
!   
   case default
      if (NEXACT.gt.0) then
         rate_error_mesh(istep) = &
         log(rel_error_mesh(istep-1)/rel_error_mesh(istep))/  &
         log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      endif
!
      rate_residual_mesh(istep) = &
      log(residual_mesh(istep-1)/residual_mesh(istep))/  &
      log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
   end select
!
!..print out the history of refinements
   write(*,*) 'HISTORY OF REFINEMENTS'
   if (NEXACT.eq.0) write(*,7005)
   if (NEXACT.gt.0) write(*,7006)
   7005  format(' mesh |','  nrdof_tot |','  nrdof_con |','    residual   |','   residual rate  ',/)
   7006  format(' mesh |','  nrdof_tot |','  nrdof_con |','    residual   |','   residual rate  |', &
            ' field error  |','rel field error|','   error rate ',/)
   !
   do i=1,istep
      if (NEXACT.eq.0) then
         write(*,7003) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_residual_mesh(i)
   7003    format(2x,i2,'  | ',2(i10,' | '),es12.5,'  |',f7.2)
      else
         write(*,7004) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_residual_mesh(i), &
                     error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
   7004    format(2x,i2,'  | ',2(i10,' | '),es12.5,'  |',f7.2,'          ', &
               2(' | ',es12.5),'  |',f7.2)
      endif
      if (i .eq. istep) write(*,*)
   enddo
!
 90 continue
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
  2010 format(A,I3,A)
   select case(Irefine)
!..uniform refinements
      case(IUNIFORM)
         if (Nreflag .eq. 1) then
!..isotropic href
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
            call global_href
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
            call update_gdof
            call update_Ddof
         elseif (Nreflag .eq. 2) then
!..anisotropic href (z)
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
            call global_href_aniso(0,1) ! refine in z
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
            call update_gdof
            call update_Ddof
         elseif (Nreflag .eq. 3) then
!..anisotropic href (z)
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
            nrelem_ref = NRELES
            do iel=1,nrelem_ref
               mdle = ELEM_ORDER(iel)
               select case(NODES(mdle)%ntype)
                  case(MDLB); kref = 1 ! refine in z
                  case(MDLP); kref = 1 ! refine in 
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
   use zoltan_wrapper
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
      call par_mumps_sc('H')
   else
      call mumps_sc('H')
   endif
!
!..do refinements and solve
   do i=0,nsteps
!..display error and refine if necessary
      if (i.ne.nsteps) then
!..Uniform refinement and solve
         call uniform_href(IUNIFORM,1,0.25d0)
         if (DISTRIBUTED .and. (.not. HOST_MESH)) then
            call par_mumps_sc('H')
         else
            call mumps_sc('H')
         endif
      else
!..Last step only display (no refinement)
         call uniform_href(INOREFINEMENT,1,0.25d0)
      endif
   enddo
!
end subroutine href_solve

