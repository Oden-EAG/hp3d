!--------------------------------------------------------------
!> Purpose : computes norm of the error and the exact solution
!!
!! @param[out]  Err   - norm of the error
!! @param[out]  Rnorm - norm of the exact solution
!--------------------------------------------------------------
   subroutine exact_error
!
      use data_structure3D
      use common_prob_data
      use environment, only : L2PROJ
!
      implicit none
!
      integer, dimension(NR_PHYSA) :: flag
!
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
 7000 format('Declare the attribute to calculate the error of: ',  &
             '   1)Displacement, 2)Stress, 3)Combined')
      read(*,*) IERROR_ATTR
!
      select case(IERROR_ATTR)
      case(DISPLACEMENT)
         flag = (/0,0,1,0/)
      case(STRESS)
         flag = (/0,0,0,1/)
      case(COMBINED)
         flag = (/0,0,1,1/)
      end select
!
      call compute_error(flag,IERROR_ATTR)
!
   end subroutine exact_error
