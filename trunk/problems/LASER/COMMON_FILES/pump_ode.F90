!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!   routine name    - pump_ode_solve
!
!----------------------------------------------------------------------
!
!   latest revision - Mar 2022
!
!   purpose         - Driver routine for computing pump power with
!                     pump ODE model, assumes pump is a plane wave
!                     over the entire core/cladding region
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
   use MPI, only: MPI_COMM_WORLD,MPI_IN_PLACE,MPI_REAL8,MPI_SUM
!
   implicit none
!
!..transverse-averaged irradiance values along the fiber
   real(8), allocatable :: pump_irr(:), sign_irr(:)
!..excited-state population density along the fiber
   real(8), allocatable :: n_ex(:)
!..auxiliary arrays
   real(8), allocatable :: zValues(:), dummy(:)
!
!..auxiliary variables
   integer :: numPts, i, j, ierr, fld
   real(8) :: a, eta, sum1, sum2, Is, Ip, dz, g0, gain
!
   if (.not. allocated(PUMP_VAL)) then
      write(*,*) 'pump_ode: PUMP_VAL has not been initiated yet. stop.'
      stop
   endif
!
!..numPts is expected to be the number of elements in z-direction
   numPts = size(PUMP_VAL)
!..dz is then the (average) element size in z-direction
   dz = ZL / numPts
!
!..allocate arrays
   allocate(pump_irr(numPts), sign_irr(numPts), n_ex(numPts))
   allocate(zValues(numPts), dummy(numPts))
!
!..fill irradiance and population density arrays
!
!  a) signal irradiance
!  a.i) compute signal power within subdomain (fiber partitioning assumed)
   a = dz/2.d0
   do i=1,numPts
      zValues(i) = (i-1)*dz+a
   enddo
!..irrationalize z values to avoid points on element interfaces
   zValues = zValues*PI*(113.d0/355.d0)
!..compute signal power in fiber core
!  (note: only sign_irr is filled with valid entries in compute_power)
   fld = 1 ! signal field index
   call compute_power(zValues,numPts,fld, dummy,pump_irr,sign_irr,n_ex)
!
!  a.ii) collect signal power values (in fiber core) on each MPI proc
   call MPI_ALLREDUCE(MPI_IN_PLACE,sign_irr,numPts,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, ierr)
!
!  a.iii) compute transverse-averaged signal irradiance over the core area
!         from signal optical power in the core area
   do i=1,numPts
      sign_irr(i) = sign_irr(i) / (PI*R_CORE*R_CORE)
   enddo
!
!  b) pump irradiance
   do i=1,numPts
      pump_irr(i) = PUMP_VAL(i) / (PI*R_CLAD*R_CLAD)
   enddo
!
!..iterate a few times (fixed-point iteration)
!  since n_ex(i) on the RHS depends on pump_irr(i)
   do j=1,5
!     c) excited-state population density
      do i=1,numPts
         Is = sign_irr(i)
         Ip = pump_irr(i)
         sum1 = (SIGMA_S_ABS/OMEGA_SIGNAL)*Is+(SIGMA_P_ABS/OMEGA_PUMP)*Ip
         sum2 = ((SIGMA_S_ABS+SIGMA_S_EMS)/OMEGA_SIGNAL)*Is + &
                ((SIGMA_P_ABS+SIGMA_P_EMS)/OMEGA_PUMP)*Ip
         n_ex(i) = sum1/(TAU_0+sum2)
      enddo
!
!  ...non-dimensional scaling factor for gain function
      g0 = ACTIVE_GAIN*L_0*SIGMA_0*NU_0
!
!  ...TODO: add option for counter-pumping, gain tailoring
!  ...solve the pump ODE by explicit stepping in z-direction
!     (pos. z-direction: co-pumped; neg. z-direction: counter-pumped)
      if (COPUMP.eq.1) then
         do i=1,numPts-1
            gain = -SIGMA_P_ABS + (SIGMA_P_ABS+SIGMA_P_EMS)*n_ex(i)
            pump_irr(i+1) = pump_irr(i) + (R_CORE*R_CORE/(R_CLAD*R_CLAD)) * &
                                        dz * g0 * gain * N_TOTAL * pump_irr(i)
         enddo
      elseif (COPUMP.eq.0) then
         do i=numPts,2,-1
            gain = -SIGMA_P_ABS + (SIGMA_P_ABS+SIGMA_P_EMS)*n_ex(i)
            pump_irr(i-1) = pump_irr(i) + (R_CORE*R_CORE/(R_CLAD*R_CLAD)) * &
                                        dz * g0 * gain * N_TOTAL * pump_irr(i)
         enddo
      else
         write(*,*) ' pump_ode_solve: COPUMP must be 1 or 0. stop.'
         stop
      endif
   enddo
!
!..update global pump power array based on irradiance solution
   do i=1,numPts
      PUMP_VAL(i) = pump_irr(i) * (PI*R_CLAD*R_CLAD)
   enddo
!
!..deallocate auxiliary arrays
   deallocate(pump_irr, sign_irr, n_ex, zValues, dummy)
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
