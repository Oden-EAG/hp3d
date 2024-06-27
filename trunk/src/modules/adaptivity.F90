!---------------------------------------------------------------------------
!> @brief Performs adaptivity based on user-provided element error estimates
!> @note: This routine is currently unused but could be overhauled to be
!!        helpful to hp3D users
!---------------------------------------------------------------------------
!
module adaptivity
!
   use environment, only : IVERBOSE_
   use error
   use data_structure3D
!
   implicit none
!
   ! save
!
   integer, parameter :: ADAPT_FILE_OUT = 35
!
   integer :: N_LIST, IPARA =0
   integer, allocatable :: NODES_LIST(:)
   real(8), allocatable :: ERROR_LIST(:)
!
   real(8) :: GREEDY
!
   real(8) :: ERR_MAX =0.d0, ERR_MIN =0.d0
!
!..ANISO_ADAP controls closure type.
!  ANISO_ADAP = 0 leads to isotropic closure
!  ANISO_ADAP = 1 leads to anisotropic closure
   integer :: ANISO_ADAP = 1
!..Parameters for anisotropic adaptation
!..0 is for greedy strat and 1 is for Doerfler
   integer, parameter :: ADAP_STRAT = 1
!..maximum number of refinements
   integer, parameter :: max_step = 200
!..ratio for the greedy or Doerfler strategy
   real(8), parameter :: Factor = 0.75
!..Factor (when multiplied with the guaranteed rate)) that decides the
!  threshold for implementing anisotropic h or p or hp
   real(8), parameter :: Factor_max_rate = 0.25
!
   contains
!
!-------------------------------------------------------------------------
!
   subroutine h_greedy_adaptivity(Element_Error)
      use data_structure3D
!
      integer :: iel, mdle, kref, istat
      real(8) :: eta
!
!  ...1 - scalar, 2,3,4 - vector directional error
      real(8), dimension(4) :: derr
      real(8), dimension(4) :: dnorm
!
      interface
         subroutine Element_Error(Mdle, Derr,Dnorm)
            integer,             intent(in)  :: Mdle
            real(8),dimension(:),intent(out) :: Derr
            real(8),dimension(:),intent(out) :: Dnorm
         end subroutine Element_Error
      end interface
!
!  ...allocation
      N_LIST = NRELES
      allocate(NODES_LIST(N_LIST), ERROR_LIST(N_LIST), STAT=istat)
      if (istat.ne.0) then
         call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
      endif

!  ...store elements
      mdle = 0
      do iel=1, N_LIST
         call nelcon(mdle, mdle)
         call Element_Error(mdle, derr, dnorm)
         NODES_LIST(iel) = mdle
         ERROR_LIST(iel) = derr(1)
      enddo
!
      call sort(NODES_LIST, ERROR_LIST, N_LIST)
      eta = GREEDY*ERROR_LIST(1)
!
      iel = 1
      do while ((ERROR_LIST(iel)>eta).and.(iel<N_LIST))
         call get_isoref(NODES_LIST(iel), kref)
         call refine(NODES_LIST(iel), kref)
         if (IVERBOSE_) then
            write(*,*) 'iel, kref = ', iel, kref
         endif
         iel = iel + 1
      enddo
!
      call close_mesh
      call update_Ddof
      call update_Gdof
!
      ERR_MAX = ERROR_LIST(1)
      ERR_MIN = ERROR_LIST(N_LIST)
!
      deallocate(NODES_LIST,ERROR_LIST, STAT=istat)
      if (istat.ne.0) then
         call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
      endif
!
   end subroutine h_greedy_adaptivity
!
end module adaptivity
