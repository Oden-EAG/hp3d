!----------------------------------------------------------------------
!
!   module name        - matrices
!
!----------------------------------------------------------------------
!
!   latest revision    - Sept 17
!
!   purpose            - store element matrices for the first layer
!                        elements
!
! Can be used for linear problems with constant material data,
! to reuse stiffness and gram matrix of a slice of elements that is
! repeated in, say the z-dimension (as in the fiber).
! Then element integration only has to be done for one slice and
! can be reuseds.
!----------------------------------------------------------------------
!
module matrices
!
   use parametersDPG
!
   implicit none
!
#include "implicit_none.h"
!
!..max # of elements in the first layer
   integer, parameter :: MAXNRFL=256
!
!..# of stored elements in the first layer
   integer :: NRFL

!
!..xy coordinates of the first vertex node
   real*8, dimension(2,3,MAXNRFL) :: XYVERT

!
!..
   VTYPE, allocatable, save :: ZFL_EE(:,:,:)
   VTYPE, allocatable, save :: ZFL_EQ(:,:,:)
   VTYPE, allocatable, save :: ZFL_QQ(:,:,:)
   VTYPE, allocatable, save :: ZFL_EF(:,:,:)
   VTYPE, allocatable, save :: ZFL_FF(:,:,:)
!
!!$OMP THREADPRIVATE (NRFL)
!!$OMP THREADPRIVATE (XYVERT)
!!$OMP THREADPRIVATE(ZFL_EE)
!!$OMP THREADPRIVATE(ZFL_EQ)
!!$OMP THREADPRIVATE(ZFL_QQ)
!!$OMP THREADPRIVATE(ZFL_EF)
!!$OMP THREADPRIVATE(ZFL_FF)
!
!..order of elements
   integer, parameter :: MYP=7
!
!..matrix dimensions
   integer, parameter :: MYE = 3*MYP*(MYP+1)**2*2
   integer, parameter :: MYF = 3*MYP*(MYP+1)**2*2
   integer, parameter :: MYQ = MYP**3*6
!
end module matrices
