!-----------------------------------------------------------------------
!
!   routine name       - mg_spcg_solve
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2018
!
!   purpose            - routine solves the system A x = b using the
!                        Conjugate Gradient solver preconditioned with 
!                        the symmetric multi-grid cycle. Everything is
!                        computed using the local matrices
!
!   arguments :
!
!                      - the solution is stored in the module pcg_data
!-----------------------------------------------------------------------
!       
   subroutine mg_pcg_solve
!
   use macro_grid_info,    only: NRELES_COARSE, DLOC, NRDOF_MACRO
   use pcg_info,           only: RHS, XSOL,LOC, TOL,  CG_ITER, MAX_ITER, &
                                 compute_global_residual, local_stiff_multiply
!
   implicit none
!
#if C_MODE
   complex*16 :: r(NRDOF_MACRO), z(NRDOF_MACRO), p(NRDOF_MACRO), Ap(NRDOF_MACRO)
   complex*16 :: zr, rz, pAp, zr_new, alpha, beta 
   complex*16, external :: zdotc
   real*8,     external :: dznrm2
#else
   real*8     :: r(NRDOF_MACRO), z(NRDOF_MACRO), p(NRDOF_MACRO), Ap(NRDOF_MACRO)
   real*8     :: zr, rz, pAp, zr_new, alpha, beta
   real*8,  external :: ddot
   real*8,  external :: dnrm2

#endif   
   integer    :: n, iter
   real*8     :: rnorm
!   
!------------------------------------------------------------------------
!
   n = NRDOF_MACRO
!..calculate initial residual
   call compute_global_residual(XSOL,r)
!
!..norm of the residual
   rnorm = dznrm2(n,r,1)
   if (rnorm .lt. TOL) return
!
!..continue otherwise. Preconditioning
   call mg_preconditioner(r,z)

!..make a copy of the correction   
   p = z
   zr = zdotc(n,z,1,r,1)
! 
   do iter = 1,MAX_ITER
!
      rz = zdotc(n,r,1,z,1)
!      
!  ...calculate Ap
      call local_stiff_multiply(p,Ap)
      pAp  = zdotc(n,p,1,Ap,1)
      alpha = rz/pAp
!
      XSOL = XSOL + alpha * p
      r = r - alpha * Ap
!
      rnorm = dznrm2(n,r,1)
      write(*,1000) iter, rnorm
 1000 format(' Iter, ||r|| =               ', i6, es13.4)      

      if (rnorm .lt. TOL ) exit

!  ...continue otherwise. Preconditioning
      call mg_preconditioner(r,z)
!
      zr_new = zdotc(n,z,1,r,1)
      beta  = zr_new/zr
      p = z + beta*p
      zr = zr_new
!  ...end of loop through iterations
   enddo
!   
   CG_ITER = iter
!
!
   end subroutine mg_pcg_solve
