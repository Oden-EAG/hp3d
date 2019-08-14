#if DEBUG_MODE

!--------------------------------------------------------
!> Purpose : contains the environmental variables for UHM
module locker
  
  use kinds
  use error
  use data_structure3D

  !------------------------------------------------------
  ! ** variables
  save
  logical, private :: data_locker = .false.

  ! relavent objects
  integer, private, target, allocatable :: &
       LOCK_MDLE(:)
  
contains
  !--------------------------------------------------------  
  logical function is_data_locked()
    is_data_locked = data_locker
  end function is_data_locked

  !> Purpose : store connectivities to array
  subroutine data_lock
    implicit none
    integer :: istat, iel, mdle

    if (is_data_locked()) return
    
    allocate(LOCK_MDLE(NRELES), stat=istat)
    if (istat.ne.SUCCESS) then
       call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
    endif

    do iel=1,NRELES
       mdle = ELEM_ORDER(iel)
       NODES(mdle)%lock = iel
       LOCK_MDLE(iel)   = mdle
    end do
    
    data_locker = .true.
  end subroutine data_lock

  !--------------------------------------------------------
  !> Purpose : deallocate workspace
  subroutine data_unlock
    implicit none
    integer :: istat

    if (.not.is_data_locked()) return

    deallocate(LOCK_MDLE, stat=istat)
    if (istat.ne.SUCCESS) then
       call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
    endif

    data_locker = .false.
  end subroutine data_unlock

  !--------------------------------------------------------
  !> Purpose : deallocate workspace
  subroutine nelcon_from_locker(Iel, Mdle)
    implicit none
    integer, intent(in)  :: Iel
    integer, intent(out) :: Mdle

    if (.not.is_data_locked()) then
       call logic_error(FAILURE,__FILE__,__LINE__)
    end if

    Mdle = LOCK_MDLE(Iel)

  end subroutine nelcon_from_locker
  
end module locker

#endif
