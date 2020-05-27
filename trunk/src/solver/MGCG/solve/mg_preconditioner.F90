!
!
#include "implicit_none.h"

   subroutine mg_preconditioner(Nrgrids,r,z,n)
!
   use mg_data_structure,only: GRID, COARSE_SOLVER, NO_CSOLVE
   use macro_grid_info,  only: NRDOF_MACRO
   use pcg_info,         only: allocate_local_vectors, deallocate_local_vectors, &
                               compute_local_residuals, update_local_residuals,  &
                               update_global_residual, XSOL, G, compute_global_residual, &
                               MAX_SM_ITER
   use parameters,       only: ZERO, ZONE
!
   implicit none
!
   integer, intent(in)  :: Nrgrids,n
   VTYPE,   intent(in)  :: r(n) 
   VTYPE,  intent(out)  :: z(n) 
   VTYPE , allocatable  :: r_aux(:),z_aux(:)

   integer  :: maxit
   integer  :: maxdepth, idepth
   integer  :: m
!
!--------------------------------------------------------   
!..for now
   maxit = MAX_SM_ITER
! 
!..total number of grids   
   maxdepth = Nrgrids - 1
!   
!..calculate updated local residuals residual

!..pre-smooth
   call block_jacobi_smoother(Nrgrids,r,n,maxit,G(nrgrids)%y1)

!..Step 3 : Compute local residuals (restrict global solution to the element)

   call compute_local_residuals(Nrgrids,XSOL+G(nrgrids)%y1,n)

!..loop through all the sub master grids and apply the smoother. Each sub-master grid is
!  one level 'down'
!
   do idepth = 1,maxdepth
!
      m = NRDOF_MACRO(nrgrids-idepth)
      allocate(r_aux(m)) ; r_aux = ZERO

      call schur_compl_restriction(nrgrids-idepth,r_aux,m)
!
      call block_jacobi_smoother(nrgrids-idepth,r_aux,m, maxit,G(nrgrids-idepth)%y1)  
!
      deallocate(r_aux)
!
!  ...Step 3 : Compute updated local residuals (restrict global sol to the element level)
      call update_local_residuals(nrgrids-idepth,G(nrgrids-idepth)%y1,m)
!
!..end of loop through all the sub master grids
   enddo   
!
   if (COARSE_SOLVER /= NO_CSOLVE) then
      call coarse_grid_solve(nrgrids-maxdepth)
   else   
      G(nrgrids-maxdepth)%y2 = ZERO
   endif   
!
!..loop through all the grids and apply the smoother. Each grid is one level 'up'
   do idepth = maxdepth,0,-1
!
      m = NRDOF_MACRO(nrgrids-idepth)
      allocate(r_aux(m)) ; r_aux = ZERO
!
      call update_global_residual(nrgrids-idepth,r_aux,m)
!
!  ...post smooth on the sub-master grid
      call block_jacobi_smoother(nrgrids-idepth,r_aux,m,maxit,G(nrgrids-idepth)%y3)  
!
      deallocate(r_aux)

!  ...exit if on the finest level
      if (idepth .eq. 0) exit
!
!  ...total correction from this level is:
      allocate(z_aux(m))
      z_aux = G(nrgrids-idepth)%y1 + G(nrgrids-idepth)%y2 + G(nrgrids-idepth)%y3
!
      call schur_compl_extension(nrgrids-idepth, z_aux,m)
      deallocate(z_aux)
!
!  ...end of loop through all the grids
   enddo
!
!..total correction to the solution
   z =  G(nrgrids)%y1 + G(nrgrids)%y2 + G(nrgrids)%y3
!
!
   end subroutine mg_preconditioner
