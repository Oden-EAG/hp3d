!--------------------------------------------------------------
!> Purpose : computes norm of the error and the exact solution
!!
!! @param[out]  Err   - norm of the error
!! @param[out]  Rnorm - norm of the exact solution
!--------------------------------------------------------------
subroutine exact_error(Err,Rnorm)
  use data_structure3D
  use common_prob_data
  implicit none
  real*8,intent(out) :: Err, Rnorm
  !
  real*8 :: derr, dnorm, err_rate
  !
  integer :: iprint, mdle, nel, nr_dof, i, n
  integer :: nr_dofH, nr_dofE, nrdofV, nrdofQ

  ! variables saved
  integer, save :: nr_dof_save
  real*8 , save :: err_save
  !------------------------------------------------------------
  iprint=0

  !  initialize
  Err = 0.d0; Rnorm = 0.d0

  ! loop over ACTIVE elements
  mdle=0

  do nel=1,NRELES
    call nelcon(mdle, mdle)
    call exact_error_element(mdle, derr,dnorm)

    Err   = Err   + derr
    Rnorm = Rnorm + dnorm

    if (iprint.eq.1) then
      write(*,7004) mdle,derr,dnorm
7004  format(' exact_error: mdle,err^2,rnorm^2 = ',i7,2x,2(e12.5,2x))
    endif
  enddo

  call find_nrdof(nr_dofH, nr_dofE, nrdofV, nrdofQ)
  nr_dof = nr_dofH

  select case(IERROR_PROB)
  case(IERROR_L2,IERROR_H1)
    Err    = sqrt(Err)
    Rnorm  = sqrt(Rnorm)
  endselect

  err_rate=0.d0
  if (IVIS.gt.0) then
     if (nr_dof.gt.nr_dof_save) then
        if (Err.gt.0.d0) then
           err_rate = dlog(err_save/Err)/log(float(nr_dof_save)/nr_dof)
        endif
     endif
  endif

  err_save = Err; nr_dof_save = nr_dof

  IVIS = IVIS + 1

  RWORK(IVIS,1) = nr_dof
  RWORK(IVIS,2) = Err
  RWORK(IVIS,3) = Rnorm
  RWORK(IVIS,4) = Err/Rnorm
  RWORK(IVIS,5) = err_rate

  write(*,*)'-- Exact Error Report -- :: IERR=', IERROR_PROB
  write(*,*)'Computed nr_dof = ', nr_dof
  write(*,*)' i, nrdof, err,rnorm,rel_err,err_rate'
  do i=1,IVIS
     write(*,9999)i,RWORK(i,1:5)
9999 format(i2,' | ',f15.1, ' | ', 4(e12.5,2x))
  enddo

end subroutine exact_error

!--------------------------------------------------------------
!> Purpose : Write error to file
!!
!! @param[in]  fp   - file to dumpout
!--------------------------------------------------------------

subroutine dumpout_error_to_file(fp)
  use common_prob_data, only : IVIS, RWORK, IERROR_PROB
  implicit none
  character(len=*), intent(in) :: fp
  integer, parameter           :: ndump = 31
  integer                      :: i
  open(unit=ndump,file=fp, &
       form='formatted',access='sequential',status='unknown')
  !
  write(ndump,*) 'Norm (IERROR_PROB) = ', IERROR_PROB
  write(ndump,*) 'Number of DOF'
  write(ndump,*) (RWORK(i,1),'  ',i=1,IVIS)
  write(ndump,*) 'Absolute Error'
  write(ndump,*) (RWORK(i,2),'  ',i=1,IVIS)
  write(ndump,*) 'Norm of Solution'
  write(ndump,*) (RWORK(i,3),'  ',i=1,IVIS)
  write(ndump,*) 'Relative Error'
  write(ndump,*) (RWORK(i,4),'  ',i=1,IVIS)
  write(ndump,*) 'Relative Error Rate'
  write(ndump,*) (RWORK(i,5),'  ',i=1,IVIS)
  !
  close(ndump)
!
  do i=1,IVIS
     write(*,9999)i,RWORK(i,1:5)
9999 format(i2,' | ',f15.1, ' | ', 4(e12.5,2x))
  enddo

end subroutine dumpout_error_to_file
