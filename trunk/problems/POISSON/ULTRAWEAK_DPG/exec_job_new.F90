

subroutine exec_job_new()

    use common_prob_data
    use data_structure3D
    use MPI           , only: MPI_COMM_WORLD
    use mpi_param     , only: RANK,ROOT,NUM_PROCS
    use par_mesh      , only: EXCHANGE_DOF,distr_mesh,DISTRIBUTED,HOST_MESH
    use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_eval
 !
    implicit none
 !
    integer :: i,ierr
    real(8) :: MPI_Wtime,start_time,end_time
    integer :: iParAttr(4) = (/0,0,0,1/)
    integer,parameter :: ref_type = 1
    EXCHANGE_DOF = .false.
    !
    if(RANK .eq. ROOT) then
        write(*,*) '=================='
        write(*,*) 'exec_job: starting'
        write(*,*) '=================='
    endif
 
    !..distribute mesh initially
    call distr_mesh
    !..set Zoltan partitioner
    if (NUM_PROCS .gt. 1) then
	 if(RANK .eq. ROOT) write(*,*) "setting up the partitioner  = ",0
	 call zoltan_w_set_lb(0)
    else
	 call zoltan_w_set_lb(0)
    endif
!..first h-refinement
    call MPI_BARRIER (MPI_COMM_WORLD, ierr);start_time = MPI_Wtime()
    ! if(RANK .eq. ROOT) 
    call global_href
    !!call close_mesh
    call MPI_BARRIER (MPI_COMM_WORLD, ierr)
    call update_gdof
    call update_Ddof
    call MPI_BARRIER (MPI_COMM_WORLD, ierr)

!..distribute mesh initially
    call distr_mesh

    call MPI_BARRIER (MPI_COMM_WORLD, ierr)
    call par_verify
    call MPI_BARRIER (MPI_COMM_WORLD, ierr)

!..second h-refinement
    !call global_href
    !!call close_mesh
    !call MPI_BARRIER (MPI_COMM_WORLD, ierr)
    !call update_gdof
    !call update_Ddof
    !call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!..distribute
    !call distr_mesh
    !call MPI_BARRIER (MPI_COMM_WORLD, ierr)
    !call par_verify
    !call MPI_BARRIER (MPI_COMM_WORLD, ierr)
 
!..set Zoltan partitioner
    !call zoltan_w_set_lb(0)

    do i=1,IMAX

        if(RANK .eq. ROOT) write(*,100) 'Beginning iteration i = ', i

        Print *, HOST_MESH, RANK
        call MPI_BARRIER (MPI_COMM_WORLD, ierr)

        if (DISTRIBUTED .and. (.not. HOST_MESH)) then
           call par_mumps_sc('H')
           !call petsc_solve('P')
        else
           call mumps_sc('H')
        endif
        call MPI_BARRIER (MPI_COMM_WORLD, ierr)

        !call exact_error
	
	!Anisotropic hp refinement
	 if(ref_type .eq. 1) then
	     call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call HpAdapt
	     call MPI_BARRIER (MPI_COMM_WORLD, ierr)
	 else
	!isotropic refinement
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call refine_DPG
	 !call uniform_href(1,1,0.75)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
	 endif
	! mesh distribution
	 if (NUM_PROCS .gt. 1) then

	    if(RANK .eq. ROOT) write(*,*) " Redistributing the mesh "

         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call distr_mesh
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
    ! mesh verification
         call par_verify
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)

     endif

    enddo

    !solving on the last mesh
    if (DISTRIBUTED .and. (.not. HOST_MESH)) then
        call par_mumps_sc('H')
        !call petsc_solve('P')
     else
        call mumps_sc('H')
     endif
    !  call exact_error

     100 format(/,'/////////////////////////////////////////////////////////////', &
     /,'             ',A,I2,/)
     call MPI_BARRIER (MPI_COMM_WORLD, ierr);end_time   = MPI_Wtime()

     iParAttr(1:4) = (/1,1,1,3/) ! write field output only
     call my_paraview_driver(iParAttr)
     call MPI_BARRIER (MPI_COMM_WORLD, ierr)
     
     if(RANK .eq. ROOT) then
        write(*,*)
        write(*,*) '=================='
        write(*,*) 'exec_job: finished'
        write(*,*) '=================='
        write(*,*)
     endif

 end subroutine exec_job_new
