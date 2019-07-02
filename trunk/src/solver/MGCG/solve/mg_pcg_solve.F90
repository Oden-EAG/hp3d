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
!
#include 'implicit_none.h'
!
   subroutine mg_pcg_solve(Igrid)
!
   use mg_data_structure,  only: GRID
   use macro_grid_info,    only: NRDOF_MACRO
   use pcg_info,           only: RHS, XSOL, TOL,  CG_ITER, MAX_ITER, &
                                 compute_global_residual, global_stiff_multiply, &
                                 allocate_local_vectors, deallocate_local_vectors

!
   implicit none
!
   integer, intent(in) :: Igrid
   VTYPE, allocatable  :: r(:), z(:), p(:), Ap(:)
   VTYPE   :: zr, rz, pAp, zr_new, alpha, beta 
   integer :: n, iter
   real*8  :: rnorm
!
   interface
      function l2norm(v)
      VTYPE,   intent(in) :: v(:)
      real*8              :: l2norm
      end function l2norm
!
      function dotp(v,u)
      VTYPE,   intent(in) :: v(:),u(:)
      VTYPE               :: dotp
      end function dotp
   end interface   
!   
!------------------------------------------------------------------------
!
   n = NRDOF_MACRO(Igrid)
   allocate(r(n),z(n),p(n),Ap(n))
!   
   call allocate_local_vectors(Igrid)
!..calculate initial residual
   call compute_global_residual(Igrid,XSOL,r,n)
!
!..norm of the residual

   rnorm = l2norm(r)

   write(*,1001)   TOL, rnorm
 1001 format(' mg_pcg_solve: tol,  ||r|| = ', es7.0, es13.4)      
   if (rnorm .lt. TOL)  then 
      write(*,1002) 
 1002 format(' mg_pcg_solve: exiting...')   
      call deallocate_local_vectors(Igrid)
      return
   endif   
!
!..continue otherwise. Preconditioning
   call mg_preconditioner(Igrid,r, z,n)

!..make a copy of the correction   
   p = z
   zr = dotp(z,r)
! 
   do iter = 1,MAX_ITER
!
      rz = dotp(r,z)
!      
!  ...calculate Ap
      call global_stiff_multiply(Igrid,p,Ap,n)
!
      pAp  = dotp(p,Ap)
      alpha = rz/pAp
      XSOL = XSOL + alpha * p
      r = r - alpha * Ap
!
      rnorm = l2norm(r)
!      
      write(*,1003) iter, rnorm
 1003 format(' mg_pcg_solve: Iter, ||r|| = ', i7, es13.4)      
!
      if (rnorm .lt. TOL ) exit
!
!  ...continue otherwise. Preconditioning
      call mg_preconditioner(Igrid,r,z,n)
!
      zr_new = dotp(z,r)
      beta  = zr_new/zr
      p = z + beta*p
      zr = zr_new
!  ...end of loop through iterations
   enddo
!   

   call deallocate_local_vectors(Igrid)

   deallocate(r,z,p,Ap)

   CG_ITER = iter
!
   end subroutine mg_pcg_solve
  
