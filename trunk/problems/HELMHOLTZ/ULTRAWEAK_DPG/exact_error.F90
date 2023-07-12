!
!----------------------------------------------------------------------
!> @brief       Compute and print exact error
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine exact_error
!  
      use data_structure3D
      use common_prob_data_UW
      use environment
      use assembly_sc,        only: NRDOF_TOT, NRDOF_CON
      use par_mesh,           only: DISTRIBUTED
      use mpi_param
      use MPI
!   
      implicit none
!
      integer :: iflag(3)
!
      real(8) :: errorH, rnormH
      real(8) :: errorE, rnormE
      real(8) :: errorV, rnormV
      real(8) :: errorQ, rnormQ
      real(8) :: err, rnorm
      integer :: iel, mdle, ierr
!
!----------------------------------------------------------------------
!
!  ...flag components to compute error for
      iflag(1) = 0; iflag(2) = 0 ; iflag(3) = 1
!
      err = 0.d0; rnorm = 0.d0
!
      if (.not. DISTRIBUTED) then
         ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
         NRELES_SUBD = NRELES
      endif
!
!$OMP PARALLEL DEFAULT(PRIVATE)        &
!$OMP SHARED(NRELES,iflag,mdle)        &
!$OMP REDUCTION(+:err,rnorm)
!$OMP DO
      do iel=1,NRELES_SUBD
         mdle = ELEM_SUBD(iel)
         call element_error(mdle,iflag,                     &
                            errorH,errorE,errorV,errorQ,    &
                            rnormH,rnormE,rnormV,rnormQ)
         err = err + errorQ
         rnorm = rnorm + rnormQ
      enddo
!$OMP END DO        
!$OMP END PARALLEL
!
      if (DISTRIBUTED) then
         call MPI_ALLREDUCE(MPI_IN_PLACE,err,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,rnorm,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif
!
      if (RANK.le.ROOT) then
         write(*,7020) NRDOF_TOT, sqrt(err),sqrt(rnorm)
         write(*,7030) NRDOF_TOT, sqrt(err/rnorm)
      endif
!
7020 format('exact_error: TOTAL NUMBER OF DOF, L2 ERROR AND NORM = ',i8,3x,2e12.5)
7030 format('exact_error: TOTAL NUMBER OF DOF, RELATIVE L2 ERROR = ',i8,3x,f8.4)
!
   end subroutine exact_error





!----------------------------------------------------------------------
!> @brief       Compute and print DPG residual
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine compute_residual
!  
      use data_structure3D
      use common_prob_data_UW
      use environment
      use assembly_sc,        only: NRDOF_TOT, NRDOF_CON
      use par_mesh,           only: DISTRIBUTED
      use mpi_param
      use MPI
!   
      implicit none
!
      real(8), allocatable :: elem_resid(:)
      integer, allocatable :: elem_ref_flag(:)
!
      real(8) :: residual
      integer :: iel, mdle, ierr
!
!----------------------------------------------------------------------
!
      allocate(elem_resid(NRELES), elem_ref_flag(NRELES))
      residual = 0.d0;
!
      if (.not. DISTRIBUTED) then
         ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
         NRELES_SUBD = NRELES
      endif
!   
!$OMP PARALLEL DEFAULT(PRIVATE)                &
!$OMP SHARED(elem_resid,elem_ref_flag,NRELES)  &
!$OMP REDUCTION(+:residual)           
!$OMP DO
      do iel=1,NRELES_SUBD
         mdle = ELEM_SUBD(iel)
         call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
         residual = residual + elem_resid(iel)
      enddo
!$OMP END DO        
!$OMP END PARALLEL
!
      if (DISTRIBUTED) then
         call MPI_ALLREDUCE(MPI_IN_PLACE,residual,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif
!
      if (RANK.le.ROOT) then
         write(*,7030) NRDOF_TOT, sqrt(residual)
      endif
 7030 format('compute_residual: TOTAL NUMBER OF DOF, RESIDUAL = ',i8,3x,1e12.5)
!
      deallocate(elem_resid, elem_ref_flag)
!
   end subroutine compute_residual
