

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
    call zoltan_w_set_lb(0)

    call MPI_BARRIER (MPI_COMM_WORLD, ierr);start_time = MPI_Wtime()
    ! if(RANK .eq. ROOT) 
    call global_href
    call update_gdof
    call update_Ddof

    call MPI_BARRIER (MPI_COMM_WORLD, ierr);end_time   = MPI_Wtime()

!..distribute mesh initially
    call distr_mesh
!..set Zoltan partitioner
    call zoltan_w_set_lb(0)

    do i=1,IMAX

        if(RANK .eq. ROOT) write(*,100) 'Beginning iteration i = ', i

        Print *, HOST_MESH, RANK
        if (DISTRIBUTED .and. (.not. HOST_MESH)) then
           call par_mumps_sc('G')
        else
           call mumps_sc('G')
        endif
        call exact_error
        call HpAdapt

    enddo

    !solving on the last mesh
    if (DISTRIBUTED .and. (.not. HOST_MESH)) then
        call par_mumps_sc('G')
     else
        call mumps_sc('G')
     endif
     call exact_error

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