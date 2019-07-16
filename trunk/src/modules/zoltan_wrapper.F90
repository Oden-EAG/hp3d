!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!     module:              zoltan_wrapper
!     last modified:       July 2019
!
!----------------------------------------------------------------------
module zoltan_wrapper
!
   use MPI
   use MPI_param
   use zoltan
!
   implicit none
!
   logical, save :: ZOLTAN_IS_INIT = .false.
!
   type(Zoltan_Struct), pointer, save :: zz
!
   contains
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_init
!     purpose:    initialize zoltan environment, and set parameters
!----------------------------------------------------------------------
   subroutine zoltan_w_init()
!
      integer(kind=Zoltan_int)   :: ierr
      real   (kind=Zoltan_float) :: ver
!
      if (ZOLTAN_IS_INIT) then
         write(*,*) 'zoltan_w_init: Zoltan has already been initialized.'
         return
      endif
!
!  ...initialize Zoltan environment
      ierr = Zoltan_Initialize(ver)
      call zoltan_w_handle_err(ierr,'Zoltan_Initialize')
!
!  ...create Zoltan memory, and set default parameters
      zz = Zoltan_Create(MPI_COMM_WORLD)
      if (.not. associated(zz)) then
         write(*,*) 'zoltan_w_init: Fatal error in Zoltan_Create.'
         return
      endif
!
!  ...set number of local ID entries to zero
      ierr = Zoltan_Set_Param(zz,'NUM_LID_ENTRIES','0')
      call zoltan_w_handle_err(ierr,'Zoltan_Set_Param')
!
      ZOLTAN_IS_INIT = .true.
!
   end subroutine zoltan_w_init
!
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_finalize
!     purpose:    close zoltan environment
!----------------------------------------------------------------------
   subroutine zoltan_w_finalize()
!
      ZOLTAN_IS_INIT = .false.
!
   end subroutine zoltan_w_finalize
!
!
!----------------------------------------------------------------------
!     routine:    zoltan_w_handle_err
!     purpose:    handle error code returned by a Zoltan function
!----------------------------------------------------------------------
   subroutine zoltan_w_handle_err(ierr,str)
      integer(kind=Zoltan_int), intent(in) :: ierr
      character(len=*)        , intent(in) :: str
      if (ierr .ne. ZOLTAN_OK) then
         write(*,*) str,': ierr = ', ierr
      endif
      ! ZOLTAN_WARN
      ! ZOLTAN_FATAL
      ! ZOLTAN_MEMERR
   end subroutine zoltan_w_handle_err
!
!
end module zoltan_wrapper
