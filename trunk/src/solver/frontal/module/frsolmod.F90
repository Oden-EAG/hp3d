!> @brief Defines space for frontal solver
!> @date Mar 2023
module frsolmod
!
   save
!..destination vectors and pointer
   integer, allocatable          :: IDESVE(:)
   integer, allocatable          :: NDESVE(:,:)
   integer                       :: NPDESV
!
!..work space for the prefront
   integer, allocatable          :: IN(:)
   integer, allocatable          :: IAWORK(:)
!
!..arrays to reorder elements in order to minimize the bandwidth
   integer, allocatable          :: NEW_ELEM_ORDER(:)
   real(8), allocatable          :: ELEM_CENTER(:)
!
!..work space for the frontal solver
   integer                       :: MFRSOL
#if C_MODE
   complex(8), allocatable       :: ZWORKFRS(:)
#else
   real(8)   , allocatable       :: ZWORKFRS(:)
#endif
!
   logical                       :: REORDER = .false.
!
contains
!
!> @brief Sets the workspace for the frontal solver
!> @param Iwork - size of workspace
   subroutine set_frsol_workspace(Iwork)
      implicit none
      integer, intent(in) :: Iwork
      MFRSOL = iwork
   end subroutine set_frsol_workspace
!
end module frsolmod
