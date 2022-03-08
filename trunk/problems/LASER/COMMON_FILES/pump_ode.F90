!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!   routine name    - pump_ode_solve
!
!----------------------------------------------------------------------
!
!   latest revision - Nov 2021
!
!   purpose         - Driver routine for computing pump power with
!                     pump ODE model
!
!   arguments
!                   - NumPts
!
!----------------------------------------------------------------------
!
subroutine pump_ode_solve
!
   use commonParam
   use laserParam
!
   implicit none
!
   integer :: numPts
!
   if (.not. allocated(PUMP_VAL)) then
      write(*,*) 'pump_ode: PUMP_VAL has not been initiated yet. stop.'
      stop
   endif
!
   numPts = size(PUMP_VAL)
!
!..Solve the pump ODE
   PUMP_VAL(1:numPts) = PLANE_PUMP_POWER
!
end subroutine pump_ode_solve
!
!----------------------------------------------------------------------
!
subroutine pump_ode_alloc(NumPts)
!
   use commonParam
   use laserParam
   use mpi_param, only: RANK,ROOT
   use MPI      , only: MPI_COMM_WORLD,MPI_IN_PLACE,MPI_REAL8,MPI_SUM
   use par_mesh , only: DISTRIBUTED,HOST_MESH
!
   implicit none
!
   integer, intent(in) :: NumPts
!
   if (allocated(PUMP_VAL)) then
      write(*,*) 'pump_ode: PUMP_VAL has already been initiated.'
   else
      allocate(PUMP_VAL(NumPts))
   endif
!
!..Initiate pump power values
   PUMP_VAL(1:NumPts) = PLANE_PUMP_POWER
!
end subroutine pump_ode_alloc
!
!----------------------------------------------------------------------
!
subroutine pump_ode_dealloc
!
   use commonParam
   use laserParam
!
   implicit none
!
   if (.not. allocated(PUMP_VAL)) then
      write(*,*) 'pump_ode_dealloc: PUMP_VAL had not been initiated.'
   else
      deallocate(PUMP_VAL)
   endif
!
end subroutine pump_ode_dealloc
!
!----------------------------------------------------------------------
