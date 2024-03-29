!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!> @name       mpi_wrapper
!> @date       Feb 2024
!----------------------------------------------------------------------
module mpi_wrapper
!
#if HP3D_USE_MPI_F08
   use mpi_f08
#else
   use MPI
#endif
   use mpi_param
   use environment, only: QUIET_MODE
!
   implicit none
!
   logical, save :: MPI_IS_INIT = .false.
!
   contains
!
!----------------------------------------------------------------------
!> @name    mpi_w_init
!> @brief   initialize MPI environment, and set parameters
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
      if (.not.QUIET_MODE .and. RANK.eq.ROOT) then
         write(*,100) 'MPI initialized successfully. NUM_PROCS = ',NUM_PROCS
      endif
  100 format(/,A,I6)
!
!  ...initialize Zoltan environment
      call zoltan_ext_init
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
!> @name    mpi_w_finalize
!> @brief   close MPI environment
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
      call zoltan_ext_finalize
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
!> @name    mpi_w_handle_err
!> @brief   handle error code returned by an MPI function
!----------------------------------------------------------------------
   subroutine mpi_w_handle_err(Ierr,Str)
      integer         , intent(in) :: Ierr
      character(len=*), intent(in) :: Str
      character(len=MPI_MAX_ERROR_STRING) :: errStr
      integer :: ierr1,len
      if (Ierr .ne. MPI_SUCCESS) then
         call MPI_Error_string(Ierr, errStr, len, ierr1)
         write(*,*) Str,': Ierr = ', Ierr, errStr(1:len)
      endif
   end subroutine mpi_w_handle_err
!
end module mpi_wrapper

!----------------------------------------------------------------------
!> @name       mpif90_wrapper
!> @brief            MPI wrapper module using F90 MPI binding
!                          (some dependencies require F90 MPI)
!> @date       Feb 2024
!----------------------------------------------------------------------
module mpif90_wrapper
!
   use mpi_param
   use MPI
   use environment, only: QUIET_MODE
!
   implicit none
!
   contains
!
!----------------------------------------------------------------------
!> @name    mpif90_w_handle_err
!> @brief   handle error code returned by an MPI function
!----------------------------------------------------------------------
   subroutine mpi_w_handle_err(Ierr,Str)
      integer         , intent(in) :: Ierr
      character(len=*), intent(in) :: Str
      character(len=MPI_MAX_ERROR_STRING) :: errStr
      integer :: ierr1,len
      if (Ierr .ne. MPI_SUCCESS) then
         call MPI_Error_string(Ierr, errStr, len, ierr1)
         write(*,*) Str,': Ierr = ', Ierr, errStr(1:len)
      endif
   end subroutine mpi_w_handle_err

end module mpif90_wrapper
