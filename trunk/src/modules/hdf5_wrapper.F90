!----------------------------------------------------------------------
!> @name       hdf5_wrapper
!> @date       Oct 2019
!----------------------------------------------------------------------
module hdf5_wrapper
!
   use HDF5
   use mpi_wrapper
!
   implicit none
!
   logical, save :: HDF5_IS_INIT = .false.
!
   contains
!
!----------------------------------------------------------------------
!> @name     hdf5_w_init
!> @brief    initialize HDF5 environment, and set parameters
!----------------------------------------------------------------------
   subroutine hdf5_w_init()
!
      integer :: ierr
!
      if (HDF5_IS_INIT) then
         write(*,*) 'HDF5_IS_INIT: HDF5 has already been initialized.'
         return
      endif
!
!  ...HDF5 initialization
      call h5open_f(ierr)
      call hdf5_w_handle_err(ierr,'h5open_f')
!
      if (RANK .eq. ROOT) then
         write(*,100) 'HDF5 initialized successfully.'
      endif
  100 format(/,A)
!
!  ...set initialization flag
      HDF5_IS_INIT = .true.
!
   end subroutine hdf5_w_init
!
!
!----------------------------------------------------------------------
!> @name    hdf5_w_finalize
!> @brief   close HDF5 environment
!----------------------------------------------------------------------
   subroutine hdf5_w_finalize()
      integer :: ierr
!
      if (.not. HDF5_IS_INIT) then
         write(*,*) 'hdf5_w_finalize: HDF5 has not been initialized.'
         return
      endif
!
!  ...close HDF5 environment
      call h5close_f (ierr)
      call hdf5_w_handle_err(ierr,'h5close_f')
!
      HDF5_IS_INIT = .false.
!
      if (RANK .eq. ROOT) then
         write(*,100) 'HDF5 closed successfully.'
      endif
  100 format(/,A)
!
   end subroutine hdf5_w_finalize
!
!----------------------------------------------------------------------
!> @name    hdf5_w_handle_err
!> @brief   handle error code returned by an HDF5 function
!----------------------------------------------------------------------
   subroutine hdf5_w_handle_err(Ierr,Str)
      integer         , intent(in) :: Ierr
      character(len=*), intent(in) :: Str
      if (Ierr .lt. 0) then
         write(*,*) str,': Ierr = ', Ierr
      endif
   end subroutine hdf5_w_handle_err
!
end module hdf5_wrapper
