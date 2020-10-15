!----------------------------------------------------------------------------
!> Purpose : save matrix to file in ascii octave format
!!
!! @param[in]  Matrix   - Matrix to save to file
!! @param[in]  Nrow     - Number of rows of matrix 
!! @param[in]  Ncol     - Number of columns of matrix 
!! @param[in]  Filename - Filename in which matrix is saved
!----------------------------------------------------------------------------
!
subroutine matrix2octave(Matrix,Nrow,Ncol,Filename)
  implicit none
  ! --------------------------------------------------------------
#if C_MODE
  complex*16,dimension(Nrow,Ncol),intent(in) :: Matrix
#else
  real*8,dimension(Nrow,Ncol),intent(in) :: Matrix
#endif
  integer,intent(in) :: Nrow,Ncol
  character(len=*), intent(in) :: Filename
  ! --------------------------------------------------------------
  integer :: nunit,ios,ii,jj
  character(len=10) :: iom
  ! --------------------------------------------------------------

  ! opening file
  nunit=100
  write(*,*) 'matrix2octave: opening file=',filename
  open(unit=nunit,iostat=ios,iomsg=iom, file=filename, &
       status="replace",action="write",form="formatted")

  ! octave reads by rows (not columns as fortran)
  write(nunit,500)
  write(nunit,600)
  write(nunit,700) 
  write(nunit,800) Nrow
  write(nunit,900) Ncol
  write(nunit,1000) ((Matrix(ii,jj), jj=1,Ncol),ii=1,Nrow)

  ! closing file
  write(*,*) 'matrix2octave: closing file=',filename
  close(unit=nunit,iostat=ios,iomsg=iom, &
       status="keep")
  
500  format('# Created by matrix2octave fortran routine')
600  format('# name: matrix')
#if C_MODE
700  format('# type: complex matrix')
#else
700  format('# type: matrix')
#endif
800  format('# rows: ',i10)
900  format('# columns: ',i10)
#if C_MODE
1000 format('(',e25.16,',',e25.16')')
#else
1000 format(e25.16)
#endif
end subroutine matrix2octave
