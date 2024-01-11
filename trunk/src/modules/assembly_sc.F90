!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!
!    module name        - assembly_sc
!
!-----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - Module sets up required workspace for the
!                         assembly of the global stiffness matrix and
!                         load vector, as well as connectivity maps to
!                         avoid recomputation after solve.
!
!    contains
!        subroutines:   - none
!
!----------------------------------------------------------------------!
module assembly_sc
!
   implicit none
!
!..timers
   real(8) :: MTime(20)
   integer, save :: IPRINT_TIME = 0
!
!..offsets for each energy space
   integer, allocatable :: NFIRSTH(:),NFIRSTE(:),NFIRSTV(:),NFIRSTQ(:)
   integer, allocatable :: NFIRST_DOF(:)
!
!..node ownership array for distributed mesh
   integer, allocatable :: NOD_OWN(:)
!
!..number of degrees of freedom
   integer :: NRDOF_CON,NRDOF_TOT
!
!..connectivity array
   integer, allocatable :: LCON(:), LCON_SUBD_CON(:)
!$OMP THREADPRIVATE (LCON,LCON_SUBD_CON)
!
!..local stiffness matrices
   VTYPE, allocatable, save :: ZLOAD(:),ZTEMP(:)
!$OMP THREADPRIVATE (ZLOAD,ZTEMP)
!
!..workspace for solution dof
   VTYPE, allocatable, save :: ZSOL_LOC(:)
!$OMP THREADPRIVATE (ZSOL_LOC)
!
end module assembly_sc
