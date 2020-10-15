!--------------------------------------------------------------
!> Purpose : Write error to file
!!
!! @param[in]  fp   - file to dumpout
!--------------------------------------------------------------

! NEED TO BE REWRITTEN

subroutine dumpout_error_to_file(fp)
!   use common_prob_data, only : IVIS, RWORK, IERROR_PROB
!   implicit none
!   character(len=*), intent(in) :: fp
!   integer, parameter           :: ndump = 31
!   integer                      :: i
!   open(unit=ndump,file=fp, &
!        form='formatted',access='sequential',status='unknown')
!   !
!   write(ndump,*) 'Norm (IERROR_PROB) = ', IERROR_PROB
!   write(ndump,*) 'Number of DOF'
!   write(ndump,*) (RWORK(i,1),'  ',i=1,IVIS)
!   write(ndump,*) 'Absolute Error'
!   write(ndump,*) (RWORK(i,2),'  ',i=1,IVIS)
!   write(ndump,*) 'Norm of Solution'
!   write(ndump,*) (RWORK(i,3),'  ',i=1,IVIS)
!   write(ndump,*) 'Relative Error'
!   write(ndump,*) (RWORK(i,4),'  ',i=1,IVIS)
!   write(ndump,*) 'Relative Error Rate'
!   write(ndump,*) (RWORK(i,5),'  ',i=1,IVIS)
!   !
!   close(ndump)
! !
!   do i=1,IVIS
!      write(*,9999)i,RWORK(i,1:5)
! 9999 format(i2,' | ',f15.1, ' | ', 4(e12.5,2x))
!   enddo

end subroutine dumpout_error_to_file