!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!     module:              mpi_wrapper
!     last modified:       July 2019
!----------------------------------------------------------------------
module mpi_wrapper
!
   use mpi_param
   use MPI, only: MPI_THREAD_FUNNELED,MPI_INIT_THREAD,         &
                  MPI_COMM_RANK,MPI_COMM_WORLD,MPI_COMM_SIZE,  &
                  MPI_BARRIER,MPI_FINALIZE,MPI_SUCCESS,        &
                  MPI_ERRORS_RETURN
   use zoltan_wrapper
!
   implicit none
!
   logical, save :: MPI_IS_INIT = .false.
!
   contains
!
!----------------------------------------------------------------------
!     routine:    mpi_w_init
!     purpose:    initialize MPI environment, and set parameters
!----------------------------------------------------------------------
   subroutine mpi_w_init()
!
      integer :: ierr,req,ret
!
      if (MPI_IS_INIT) then
         write(*,*) 'mpi_w_init: MPI has already been initialized.'
         return
      endif
!
!  ...MPI initialization
      req = MPI_THREAD_FUNNELED
      call MPI_INIT_THREAD (req, ret,ierr)
      call mpi_w_handle_err(ierr,'MPI_INIT_THREAD')
      if (ret < req) then
         write(*,*) 'mpi_wrapper: MPI_INIT_THREAD, ', &
            'provided thread level < requested thread level. stop.'
         stop
      endif
!
!  ...find out my RANK (process id), and number of processors
      call MPI_COMM_RANK (MPI_COMM_WORLD, RANK, ierr)
      call mpi_w_handle_err(ierr,'MPI_COMM_RANK')
      call MPI_COMM_SIZE (MPI_COMM_WORLD, NUM_PROCS, ierr)
      call mpi_w_handle_err(ierr,'MPI_COMM_SIZE')
!
      if (.not. QUIET_MODE .and. RANK .eq. ROOT) then
         write(*,100) 'MPI initialized successfully. NUM_PROCS = ',NUM_PROCS
      endif
  100 format(/,A,I4)
!
!  ...initialize Zoltan environment
      call zoltan_w_init
!
!  ...set initialization flag
      MPI_IS_INIT = .true.
!
!      call MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr);
!      call mpi_w_handle_err(ierr,'MPI_ERR_HANDLER')
!
   end subroutine mpi_w_init
!
!
!----------------------------------------------------------------------
!     routine:    mpi_w_finalize
!     purpose:    close MPI environment
!----------------------------------------------------------------------
   subroutine mpi_w_finalize()
      integer :: ierr
!
      if (.not. MPI_IS_INIT) then
         write(*,*) 'mpi_w_finalize: MPI has not been initialized.'
         return
      endif
!
!  ...close Zoltan environment
      call zoltan_w_finalize
!
!  ...close MPI environment
      call MPI_FINALIZE (ierr)
      call mpi_w_handle_err(ierr,'MPI_FINALIZE')
!
      MPI_IS_INIT = .false.
!
   end subroutine mpi_w_finalize
!
!----------------------------------------------------------------------
!     routine:    mpi_w_handle_err
!     purpose:    handle error code returned by an MPI function
!----------------------------------------------------------------------
   subroutine mpi_w_handle_err(Ierr,Str)
      integer         , intent(in) :: Ierr
      character(len=*), intent(in) :: Str
      character(len=64) :: errStr
      integer :: ierr1,len
      if (Ierr .ne. MPI_SUCCESS) then
         call MPI_Error_string(Ierr, errStr, len, ierr1)
         write(*,*) Str,': Ierr = ', Ierr, errStr(1:len)
      endif
   end subroutine mpi_w_handle_err
!
end module mpi_wrapper
