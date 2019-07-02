!> Purpose : module defines space for frontal solver
!! @revision Dec 10
module frsolmod
  !
  save
  !  ...destination vectors and pointer
  integer, allocatable          :: IDESVE(:)
  integer, allocatable          :: NDESVE(:,:)
  integer                       :: NPDESV
  !
  !  ...work space for the prefront
  integer, allocatable          :: IN(:)
  integer, allocatable          :: IAWORK(:)
  !
  !  ...arrays to reorder elements in order to minimize the bandwidth
  integer, allocatable          :: OLD_ELEM_ORDER(:)
  integer, allocatable          :: NEW_ELEM_ORDER(:)
  double precision, allocatable :: ELEM_CENTER(:)
  !
  !  ...work space for the frontal solver
  integer                       :: MFRSOL
#if C_MODE
  complex*16, allocatable       :: ZWORKFRS(:)
#else
  double precision, allocatable :: ZWORKFRS(:)
#endif
  !
  logical                       :: REORDER = .false.
  !
contains

  !> Purpose : set the workspace for the frontal solver
  !! @param iwork size of workspace
  subroutine set_frsol_workspace(iwork)
    implicit none
    integer, intent(in) :: iwork
    MFRSOL = iwork
  end subroutine set_frsol_workspace

end module frsolmod
