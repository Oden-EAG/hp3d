!--------------------------------------------------------------
!> Purpose : computes norm of the error and the exact solution
!!
!! @param[out]  Err   - norm of the error
!! @param[out]  Rnorm - norm of the exact solution
!--------------------------------------------------------------
subroutine exact_error
  use physics
  use common_prob_data
  use environment, only : L2PROJ
  implicit none
!
  integer, dimension(NR_PHYSA) :: flag
!------------------------------------------------------------
!
  select case(IERROR_PROB)
  case(IERROR_L2)
    L2PROJ = .TRUE.
  case(IERROR_NATURAL)
    L2PROJ = .FALSE.
  case default
    write(*,*) 'exact_error : Error calculation type not supported'
  end select
!
 333  write(*,7000)
 7000 format('Declare the attribute to calculate the error of: ')
      write(*,*)'   1)Field Displ, 2)Field Stress, 3)Field Combined, 4)Field Lagrangian, 5) Displ trace, 6) traction, 7) combined trace'
  read(*,*) IERROR_ATTR
!
  select case(IERROR_ATTR)
    case(DISPLACEMENT)
      flag = (/0,0,1,0,0/)
    case(STRESS)
      flag = (/0,0,0,1,0/)
    case(COMBINED)
      flag = (/0,0,1,1,0/)
    case(LAGRANGE)
      flag = (/0,0,0,0,1/)
    case(5)
      flag = (/1,0,0,0,0/)
    case(6)
      flag = (/0,1,0,0,0/)
    case(7)
      flag = (/1,1,0,0,0/)
  end select
!
  call compute_error(flag,IERROR_ATTR)
!
end subroutine exact_error
