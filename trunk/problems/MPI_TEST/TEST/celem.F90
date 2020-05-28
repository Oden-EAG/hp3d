!--------------------------------------------------------------------
!
!     routine name      - celem
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          - calculate modified stiffness matrix and load vector
!
!     arguments:
!
!     in:
!             Mdle      - middle node of an element
!             Idec      - 1 only info for nodes, 2 matrices
!     out:
!             Nrdofs    - # of element local dof for each physics
!             Nrdofm    - # of modified element dof in the expanded
!             Nrdofc    - # of modified element dof after compression
!             Nodm      - actual nodes in the order
!             NdofmH    - the number of H1 dof
!             NdofmE    - the number of H(curl) dof
!             NdofmV    - the number of H(div) dof
!             NdofmQ    - the number of L2 dof
!             Nrnodm    - number of the modified element nodes
!             Bload     - 1D array containing the modified load vector
!             Astif     - 1D array containing the modified stiffness matrix!
!---------------------------------------------------------------------
#include "typedefs.h"
!
   subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,Nodm, &
                    NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
   use physics
   use data_structure3D
!
!--------------------------------------------------------------------------
!
  implicit none
!
   integer,                      intent(in)  :: Mdle,Idec
   integer, dimension(NR_PHYSA), intent(out) :: Nrdofs
   integer,                      intent(out) :: Nrdofm,Nrdofc
   integer, dimension(MAXNODM),  intent(out) :: Nodm
   integer, dimension(MAXNODM),  intent(out) :: NdofmH,NdofmE,NdofmV,NdofmQ
   integer,                      intent(out) :: Nrnodm
   VTYPE  ,                      intent(out) :: Bload(*),Astif(*)
   integer, dimension(NR_PHYSA)              :: nbcond
!
!--------------------------------------------------------------------------
!
!..This is a hack to eliminate the bubbles in the trace physics variables,
!  so that only the boundary dof are included for the trace variables.
!  This is done by saying that the value is known (to be zero) for the
!  bubble (middle/interior) dof by using the BC flag (1 if known, 0 if unknown)
!  The 'known' dof are eliminated by static condensation locally
   nbcond = (/1,0/) ! removes the bubbles
   call encod(nbcond,10,NR_PHYSA, NODES(Mdle)%bcond)
   call set_index(NODES(Mdle)%case,NODES(Mdle)%bcond, NODES(Mdle)%index)
!
!..redirect to the system routine
   call celem_system(Mdle,Idec, &
                     Nrdofs,Nrdofm,Nrdofc, &
                     Nodm,NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm, &
                     Bload,Astif)
!
!..reset the BC flags back to zero (unknown)
   NODES(Mdle)%bcond = 0
   call set_index(NODES(Mdle)%case,NODES(Mdle)%bcond, NODES(Mdle)%index)
!
!
   end subroutine celem
