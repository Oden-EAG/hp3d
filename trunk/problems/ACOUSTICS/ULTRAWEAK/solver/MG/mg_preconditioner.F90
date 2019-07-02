





   subroutine mg_preconditioner(r,z)
!
   use macro_grid_info,  only: NRDOF_MACRO, NRELES_COARSE
   use pcg_info,         only: allocate_local_vectors, deallocate_local_vectors, &
                               compute_local_residuals, compute_global_residual, XSOL
   use common_prob_data                               
   use parameters,       only: ZERO
!
   implicit none
!
#if C_MODE
   complex*16, intent(in)  :: r(NRDOF_MACRO) 
   complex*16, intent(out) :: z(NRDOF_MACRO) 
   complex*16              :: xnew(NRDOF_MACRO) 
   complex*16              :: rnew(NRDOF_MACRO) 
   complex*16              :: z1(NRDOF_MACRO), z2(NRDOF_MACRO), z3(NRDOF_MACRO)   
#else
   real*8,     intent(in)  :: r(NRDOF_MACRO) 
   real*8,     intent(out) :: z(NRDOF_MACRO) 
   real*8                  :: xnew(NRDOF_MACRO) 
   real*8                  :: rnew(NRDOF_MACRO) 
   complex*16              :: z1(NRDOF_MACRO), z2(NRDOF_MACRO), z3(NRDOF_MACRO)   
#endif 
   integer  :: maxit, iel
   real*8   :: start, finish, omp_get_wtime
!
!--------------------------------------------------------   
!..for now
   maxit = 10
!   
!..pre-smooth

   ! start = omp_get_wtime()

   call block_jacobi_smoother(r,NRDOF_MACRO,maxit, z1)

   ! finish = omp_get_wtime()

   ! write(*,1000) finish-start
 ! 1000 format('presmooth time  = ',es13.4)

!   
!..calculate updated local residuals residual
   call allocate_local_vectors

!..updated solution after correction
   xnew = XSOL + z1 
!


   ! start = omp_get_wtime()
   call compute_local_residuals(xnew)
   ! finish = omp_get_wtime()

   ! write(*,1001) finish-start
 ! 1001 format('local_res time  = ',es13.4)

   if (ISOLVER == ISOLVER_MG) then
      call coarse_grid_solve(z2)
   else   
      z2 = ZERO
   endif   
!
!..update correction
   z = XSOL + z1 + z2
!
   ! start = omp_get_wtime()
   call compute_global_residual(z,rnew)
   ! finish = omp_get_wtime()
   ! write(*,1002) finish-start
 ! 1002 format('glocal_res time = ',es13.4)   
!
!..post-smooth

   
   ! start = omp_get_wtime()
   call block_jacobi_smoother(rnew,NRDOF_MACRO,maxit, z3)   
   ! finish = omp_get_wtime()
   ! write(*,1003) finish-start
 ! 1003 format('postmooth time  = ',es13.4)   
!
   z = z1 + z2 + z3
!
   call deallocate_local_vectors
!
!
   end subroutine mg_preconditioner
