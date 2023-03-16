!
#include "typedefs.h"
#include <petsc/finclude/petscksp.h>
! -----------------------------------------------------------------------
!
!    module name        - petsc_w_ksp
!
! -----------------------------------------------------------------------
!
!    latest revision    - Feb 2020
!
!    purpose            - module sets up required workspace for
!                         interfacing with linear PETSc solvers
!
! -----------------------------------------------------------------------
module petsc_w_ksp
!
   use MPI      , only: MPI_COMM_WORLD,MPI_COMM_SELF
   use mpi_param, only: RANK,NUM_PROCS
   use petscksp
!
   implicit none
!
   KSP     petsc_ksp
   KSPType petsc_ksp_type
   PC      petsc_pc
   PCType  petsc_pc_type
   Mat     petsc_A
   Vec     petsc_rhs
   Vec     petsc_sol
!
   contains
!
! -----------------------------------------------------------------------
!
subroutine petsc_ksp_start
  PetscErrorCode ierr
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
  call KSPCreate(MPI_COMM_WORLD, petsc_ksp,ierr); CHKERRQ(ierr)
  petsc_pc_type=""
end subroutine petsc_ksp_start
!
! -----------------------------------------------------------------------
!
subroutine petsc_ksp_destroy
   PetscErrorCode ierr
   call KSPDestroy(petsc_ksp, ierr); CHKERRQ(ierr)
   call PetscFinalize(ierr); CHKERRQ(ierr)
end subroutine petsc_ksp_destroy
!
end module petsc_w_ksp
