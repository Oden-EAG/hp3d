!----------------------------------------------------------------------
! exec_job_adap_ref : for adaptive refinements (May 2023)
!----------------------------------------------------------------------
subroutine exec_job_adap_ref()
!
   use common_prob_data
   use data_structure3D
   use mpi_wrapper
   use par_mesh      , only: EXCHANGE_DOF,distr_mesh,DISTRIBUTED,HOST_MESH
   use paraview      , only: paraview_select_attr
   use zoltan_wrapper, only: zoltan_w_set_lb,zoltan_w_eval
   use adaptivity    , only: ANISO_ADAP
!
   implicit none
!
   integer :: i,ierr
   real(8) :: start_time,end_time
   logical :: iPvAttr(4)
!
!----------------------------------------------------------------------
!
   EXCHANGE_DOF = .false.
!
   if(RANK .eq. ROOT) then
      write(*,*) '==========================='
      write(*,*) 'exec_job_adap_ref: starting'
      write(*,*) '==========================='
   endif
!
!..distribute mesh initially
   call distr_mesh
!
!..set Zoltan partitioner
   if (NUM_PROCS .gt. 1) then
      if(RANK .eq. ROOT) write(*,*) "setting up the partitioner"
      call zoltan_w_set_lb(0)
   endif
!
!..first h-refinement
   call MPI_BARRIER (MPI_COMM_WORLD, ierr);start_time = MPI_Wtime()
!
   call global_href
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   call update_gdof
   call update_Ddof
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..distribute mesh initially
   call distr_mesh
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   call par_verify
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
   do i=1,IMAX
!
      if(RANK .eq. ROOT) write(*,100) 'Beginning iteration i = ', i
!
      write(*,*) HOST_MESH, RANK
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
      if (DISTRIBUTED .and. (.not. HOST_MESH)) then
         call par_mumps_sc('H')
      else
         call mumps_sc('H')
      endif
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
      if(ANISO_ADAP .eq. 1) then
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call HpAdapt     
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
     else  
!..isotropic refinement
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call refine_DPG
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
     endif
!
!  ...mesh distribution
      if (NUM_PROCS .gt. 1) then
!
         if(RANK .eq. ROOT) write(*,*) " Redistributing the mesh "
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         call distr_mesh
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!     ...mesh verification
         call par_verify
         call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
      endif
   enddo
!
!..solving on the last mesh
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      call par_mumps_sc('H')
   else
      call mumps_sc('H')
   endif
!
   100 format(/,'/////////////////////////////////////////////////////////////', &
              /,'             ',A,I2,/)
   call MPI_BARRIER (MPI_COMM_WORLD, ierr);end_time   = MPI_Wtime()
!
   iPvAttr(1:4) = (/.false.,.false.,.true.,.true./) ! write field output only
   call paraview_select_attr(iPvAttr)
   call paraview_driver
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
   if(RANK .eq. ROOT) then
      write(*,*)
      write(*,*) '==========================='
      write(*,*) 'exec_job_adap_ref: finished'
      write(*,*) '==========================='
      write(*,*)
   endif
!
end subroutine exec_job_adap_ref
