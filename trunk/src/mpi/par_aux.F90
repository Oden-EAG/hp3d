!
!..auxiliary subroutines
!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!     subroutine:          collect_dofs
!
!     last modified:       July 2019
!
!     purpose:             collect solution degrees of freedom from all
!                          processors into the ROOT processor
!
!----------------------------------------------------------------------
subroutine collect_dofs()
!
   use data_structure3D
   use par_mesh , only: DISTRIBUTED
   use MPI_param, only: ROOT,RANK
   use MPI      , only: MPI_INTEGER,MPI_COMM_WORLD
!
   implicit none
!
!..auxiliary variables
   integer :: ierr
!
!----------------------------------------------------------------------
!
   if (.not. DISTRIBUTED) then
      write(*,*) 'collect_dofs: mesh is not distributed.'
      goto 190
   endif
!
!..collect degrees of freedom from every element
!..could use visit flags to avoid re-sending information
!
  190 continue
!
end subroutine collect_dofs
