!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!   routine name       - pump_ode
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
subroutine pump_ode(NumPts)
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
   if (.not. allocated(PUMP_VAL)) then
      write(*,*) 'pump_ode: PUMP_VAL has not been initiated. stop.'
      stop
   endif
!
!..Solve the pump ODE
   PUMP_VAL(1:NumPts) = PLANE_PUMP_POWER
!
end subroutine pump_ode
