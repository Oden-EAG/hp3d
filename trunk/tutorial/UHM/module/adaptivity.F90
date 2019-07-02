!---------------------------------------------------------------------------
!> Purpose : adaptivity
!---------------------------------------------------------------------------
!
module adaptivity_local
  !
  use environment, only : IVERBOSE_
  use error
  use data_structure3D
  save

  integer, parameter :: ADAPT_FILE_OUT = 35

  integer :: N_LIST, IPARA =0
  integer, allocatable :: NODES_LIST(:)
  real*8 , allocatable :: ERROR_LIST(:)

  real*8 :: GREEDY

  real*8 :: ERR_MAX =0.d0, ERR_MIN =0.d0
  !

contains
  !-------------------------------------------------------------------------
  subroutine h_adaptivity
    use data_structure3D
    implicit none
    integer :: iel, mdle, kref, istat
    real*8 :: eta

    real*8,dimension(MAXEQNH) :: derr
    real*8,dimension(MAXEQNH) :: dnorm

!     interface
!        subroutine Element_Error(Mdle, Derr, Dnorm)
!          integer,            intent(in)  :: Mdle
!          real*8,dimension(:),intent(out) :: Derr
!          real*8,dimension(:),intent(out) :: Dnorm
!        end subroutine Element_Error
!     end interface

    ! allocation
    N_LIST = NRELES
    allocate(NODES_LIST(N_LIST), ERROR_LIST(N_LIST), STAT=istat)
    if (istat.ne.0) then
       call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
    endif

    ! store elements
    mdle = 0
    do iel=1, N_LIST
       call nelcon(mdle, mdle)
       !call exact_error_element(mdle, derr, dnorm)
       call exact_error_element(mdle, derr, dnorm)
       NODES_LIST(iel) = mdle
       ERROR_LIST(iel) = derr(1)
    enddo

    call sort(NODES_LIST, ERROR_LIST, N_LIST)
    eta = GREEDY*ERROR_LIST(1)

    
    iel = 1
    do while ((ERROR_LIST(iel)>eta).and.(iel<N_LIST)) 
       call get_isoref(NODES_LIST(iel), kref)
       call refine(NODES_LIST(iel), kref)
       if (IVERBOSE_) then
          write(*,*) 'iel, kref = ', iel, kref
       end if
       iel = iel + 1
    end do
    call close
    !call update_ddof
    !call update_gdof

    ERR_MAX = ERROR_LIST(1)
    ERR_MIN = ERROR_LIST(N_LIST)

    deallocate(NODES_LIST, ERROR_LIST, STAT=istat)
    if (istat.ne.0) then
       call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
    endif

  end subroutine h_adaptivity
end module adaptivity_local
