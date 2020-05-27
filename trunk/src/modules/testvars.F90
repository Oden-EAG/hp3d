!----------------------------------------------------------------------------
!> Purpose : module defines test variables info
!!
!!   @date Jun 15
!----------------------------------------------------------------------------
module testvars
  !
  use environment , only : QUIET_MODE
  !
  implicit none
  ! the file with the test variables data
  character(len=128), save :: FILE_TESTVARS
  ! number of different test attributes
  integer, save :: NR_TEST
  ! list of test attributes
  character(len=5), save, allocatable :: TEST_NAMES(:)
  ! the corresponding number of components of each test  variable
  integer, save, allocatable :: NR_TESTCOMP(:)
  ! the corresponding discretization type
  character(len=6), save, allocatable :: TEST_SPACES(:)
  !
  interface dumpout_testvars
     module procedure dumpout_testvars_to_file
     module procedure dumpout_testvars_to_default
  end interface
  !
  interface dumpin_testvars
     module procedure dumpin_testvars_from_file
     module procedure dumpin_testvars_from_default
  end interface
  !
contains
  !
  !> Purpose : allocate data structure for test variables
  subroutine alloc_testvars
    !
    if (allocated(TEST_NAMES)) then
       deallocate(TEST_NAMES,NR_TESTCOMP,TEST_SPACES)
    endif
    !
    allocate(TEST_NAMES(NR_TEST))
    allocate(NR_TESTCOMP(NR_TEST))
    allocate(TEST_SPACES(NR_TEST))
    !
  end subroutine alloc_testvars

  !> Purpose : deallocate data structure for multiphysics
  subroutine dealloc_testvars
    if (allocated(TEST_NAMES)) then
       deallocate(TEST_NAMES,NR_TESTCOMP,TEST_SPACES)
    endif
  end subroutine dealloc_testvars
  !
  !> Purpose : routine dumps out multiphysics data structure
  !! @param fp file to dumpout
  subroutine dumpout_testvars_to_file(fp)
    implicit none
    character(len=*), intent(in) :: fp
    integer, parameter           :: ndump = 31
    integer                      :: i
    open(unit=ndump,file=fp, &
         form='formatted',access='sequential',status='unknown')
    !
    write(ndump,*) NR_TEST
    write(ndump,*) (TEST_NAMES(i),'  ',i=1,NR_TEST)
    write(ndump,*) NR_TESTCOMP
    write(ndump,*) (TEST_SPACES(i),'  ',i=1,NR_TEST)
    !
    close(ndump)
  end subroutine dumpout_testvars_to_file
  !
  !> Purpose : dumpout to default location "files/dumpTEST".
  subroutine dumpout_testvars_to_default
    call dumpout_testvars_to_file('files/dumpTEST')
  end subroutine dumpout_testvars_to_default
  !
  !> Purpose : routine dumps in test variables data structure
  !! @param fp file to dumpout
  subroutine dumpin_testvars_from_file(fp)
    implicit none
    character(len=*), intent(in) :: fp
    integer, parameter           :: ndump = 31
    integer                      :: i
    open(unit=ndump,file=fp, &
         form='formatted',access='sequential',status='unknown')
    !
    read(ndump,*) NR_TEST
    allocate(TEST_NAMES(NR_TEST))
    read(ndump,*) TEST_NAMES
    allocate(NR_TESTCOMP(NR_TEST))
    read(ndump,*) NR_TESTCOMP
    allocate(TEST_SPACES(NR_TEST))
    read(ndump,*) TEST_SPACES
    !
    close(ndump)
  end subroutine dumpin_testvars_from_file
  !
  !> Purpose : dumpin to default location "files/dumpTEST".
  subroutine dumpin_testvars_from_default
    call dumpin_testvars_from_file('files/dumpTEST')
  end subroutine dumpin_testvars_from_default
  !
  !> Purpose : read input file and set physics and data structure
  subroutine read_testvars(Fp)
    implicit none
    character(len=*), intent(in) :: Fp
    integer, parameter :: nin = 103
    integer :: i
    logical :: exist
    !----------------------------------------------------------------------
    ! file
    inquire(file=Fp,exist=exist)
    if (.NOT.exist) then
      write(*,*) ''
      write(*,*) 'File for test variables does not exist.'
      write(*,*) 'Proceeding without this information.'
      write(*,*) ''
    else
      open(unit=nin,file=Fp, &
            form='formatted',access='sequential',status='unknown')
      !  ...read number of test variables
      read(nin,*) NR_TEST
      !  ...allocate the data structure arrays
      call alloc_testvars
      !  ...read in the test attributes, the corresponding discretization
      !     type and number of components
      do i=1,NR_TEST
        read(nin,*) TEST_NAMES(i),TEST_SPACES(i),NR_TESTCOMP(i)
      enddo
      close(nin)
    !----------------------------------------------------------------------
      ! printing
      if (.NOT. QUIET_MODE) write(*,*)''
      if (.NOT. QUIET_MODE) write(*,*) '-- Test Variables --'
      if (.NOT. QUIET_MODE) write(*,*) 'TEST_NAMES | TEST_SPACES | NR_TESTCOMP'
      do i=1, NR_TEST
        if (.NOT. QUIET_MODE) then
          write(*,1003) TEST_NAMES(i),TEST_SPACES(i),NR_TESTCOMP(i)
1003 format(3x,a5,9x,a6,10x,i1)
        endif
      enddo
      if (.NOT. QUIET_MODE) write(*,*)''
    endif
    !
  end subroutine read_testvars
!
end module testvars
