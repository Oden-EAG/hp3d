
module matrix2octave_module

  implicit none
  save

  interface matrix2octave
     module procedure matrix2octave_real
     module procedure vector2octave_real
     module procedure matrix2octave_complex
     module procedure vector2octave_complex
  end interface

contains


  !----------------------------------------------------------------------------
  !> Purpose : save matrix to file in ascii octave format
  !!
  !! @param[in]  Matrix   - Matrix to save to file
  !! @param[in]  Nrow     - Number of rows of matrix
  !! @param[in]  Ncol     - Number of columns of matrix
  !! @param[in]  Filename - Filename in which matrix is saved
  !----------------------------------------------------------------------------
  !
  subroutine matrix2octave_real(Matrix,Nrow,Ncol,Filename)
    implicit none
    ! --------------------------------------------------------------
    real(kind=kind(1.d0)),dimension(Nrow,Ncol),intent(in) :: Matrix
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

500 format('# Created by matrix2octave fortran routine')
600 format('# name: matrix')
700 format('# type: matrix')
800 format('# rows: ',i10)
900 format('# columns: ',i10)
1000 format(e25.16)

  end subroutine matrix2octave_real

  !----------------------------------------------------------------------------
  !> Purpose : save vector to file in ascii octave format
  !!
  !! @param[in]  Vector   - Vector to save to file
  !! @param[in]  Nrow     - Number of rows of vector
  !! @param[in]  Filename - Filename in which vector is saved
  !----------------------------------------------------------------------------
  !
  subroutine vector2octave_real(Vector,Nrow,Filename)
    implicit none
    ! --------------------------------------------------------------
    real(kind=kind(1.d0)),dimension(Nrow),intent(in) :: Vector
    integer,intent(in) :: Nrow
    character(len=*), intent(in) :: Filename
    ! --------------------------------------------------------------
    integer :: ncol
    integer :: nunit,ios,ii
    character(len=10) :: iom
    ! --------------------------------------------------------------
    ncol=1
    ! opening file
    nunit=100
    write(*,*) 'matrix2octave: opening file=',filename
    open(unit=nunit,iostat=ios,iomsg=iom, file=filename, &
         status="replace",action="write",form="formatted")

    ! writing into file
    write(nunit,500)
    write(nunit,600)
    write(nunit,700)
    write(nunit,800) Nrow
    write(nunit,900) ncol
    write(nunit,1000) (Vector(ii), ii=1,Nrow)

    ! closing file
    write(*,*) 'matrix2octave: closing file=',filename
    close(unit=nunit,iostat=ios,iomsg=iom, &
         status="keep")

500 format('# Created by matrix2octave fortran routine')
600 format('# name: matrix')
700 format('# type: matrix')
800 format('# rows: ',i10)
900 format('# columns: ',i10)
1000 format(e25.16)

  end subroutine vector2octave_real

  !----------------------------------------------------------------------------
  !> Purpose : save matrix to file in ascii octave format
  !!
  !! @param[in]  Matrix   - Matrix to save to file
  !! @param[in]  Nrow     - Number of rows of matrix
  !! @param[in]  Ncol     - Number of columns of matrix
  !! @param[in]  Filename - Filename in which matrix is saved
  !----------------------------------------------------------------------------
  !
  subroutine matrix2octave_complex(Matrix,Nrow,Ncol,Filename)
    implicit none
    ! --------------------------------------------------------------
    complex(kind=kind(1.d0)),dimension(Nrow,Ncol),intent(in) :: Matrix
    integer,intent(in) :: Nrow,Ncol
    character(len=*), intent(in) :: Filename
    ! --------------------------------------------------------------
    integer :: nunit_real,nunit_imag,ios,ii,jj
    character(len=10) :: iom
    character(len=len(Filename)+5) :: filename_real,filename_imag
    ! --------------------------------------------------------------

    ! Don't know how to save complex number in octave ascii format ->
    ! saving real part and imaginary parts in separate files

    ! opening file real part
    nunit_real=100
    filename_real='real_'//Filename
    write(*,*) 'matrix2octave: opening files=',filename_real
    open(unit=nunit_real,iostat=ios,iomsg=iom, file=filename_real, &
         status="replace",action="write",form="formatted")

    ! opening file imag part
    nunit_imag=102
    filename_imag='imag_'//Filename
    write(*,*) 'matrix2octave: opening files=',filename_imag
    open(unit=nunit_imag,iostat=ios,iomsg=iom, file=filename_imag, &
         status="replace",action="write",form="formatted")

    ! Note that octave reads by rows (not columns as fortran)

    ! writing the real part to file
    write(nunit_real,500)
    write(nunit_real,600)
    write(nunit_real,700)
    write(nunit_real,800) Nrow
    write(nunit_real,900) Ncol
    write(nunit_real,1000) ((real(Matrix(ii,jj)), jj=1,Ncol),ii=1,Nrow)

    ! writing the imag part to file
    write(nunit_imag,500)
    write(nunit_imag,600)
    write(nunit_imag,700)
    write(nunit_imag,800) Nrow
    write(nunit_imag,900) Ncol
    write(nunit_imag,1000) ((aimag(Matrix(ii,jj)), jj=1,Ncol),ii=1,Nrow)

    ! closing files
    write(*,*) 'matrix2octave: closing file=',filename_real
    close(unit=nunit_real,iostat=ios,iomsg=iom, &
         status="keep")

    write(*,*) 'matrix2octave: closing file=',filename_imag
    close(unit=nunit_imag,iostat=ios,iomsg=iom, &
         status="keep")


500 format('# Created by matrix2octave fortran routine')
600 format('# name: matrix')
700 format('# type: matrix')
800 format('# rows: ',i10)
900 format('# columns: ',i10)
1000 format(e25.16)
    !
!!$1000 format('(',e25.16,',',e25.16')')

  end subroutine matrix2octave_complex

  !----------------------------------------------------------------------------
  !> Purpose : save vector to file in ascii octave format
  !!
  !! @param[in]  Vector   - Vector to save to file
  !! @param[in]  Nrow     - Number of rows of vector
  !! @param[in]  Filename - Filename in which vector is saved
  !----------------------------------------------------------------------------
  !
  subroutine vector2octave_complex(Vector,Nrow,Filename)
    implicit none
    ! --------------------------------------------------------------
    complex(kind=kind(1.d0)),dimension(Nrow),intent(in) :: Vector
    integer,intent(in) :: Nrow
    character(len=*), intent(in) :: Filename
    ! --------------------------------------------------------------
    integer :: ncol
    integer :: nunit_real,nunit_imag,ios,ii
    character(len=10) :: iom
    character(len=len(Filename)+5) :: filename_real,filename_imag
    ! --------------------------------------------------------------
    ncol=1
    ! Don't know how to save complex number in octave ascii format ->
    ! saving real part and imaginary parts in separate files

    ! opening file real part
    nunit_real=100
    filename_real='real_'//Filename
    write(*,*) 'matrix2octave: opening files=',filename_real
    open(unit=nunit_real,iostat=ios,iomsg=iom, file=filename_real, &
         status="replace",action="write",form="formatted")

    ! opening file imag part
    nunit_imag=102
    filename_imag='imag_'//Filename
    write(*,*) 'matrix2octave: opening files=',filename_imag
    open(unit=nunit_imag,iostat=ios,iomsg=iom, file=filename_imag, &
         status="replace",action="write",form="formatted")

    ! writing the real part to file
    write(nunit_real,500)
    write(nunit_real,600)
    write(nunit_real,700)
    write(nunit_real,800) Nrow
    write(nunit_real,900) ncol
    write(nunit_real,1000) (real(Vector(ii)),ii=1,Nrow)

    ! writing the imag part to file
    write(nunit_imag,500)
    write(nunit_imag,600)
    write(nunit_imag,700)
    write(nunit_imag,800) Nrow
    write(nunit_imag,900) ncol
    write(nunit_imag,1000) (aimag(Vector(ii)),ii=1,Nrow)

    ! closing files
    write(*,*) 'matrix2octave: closing file=',filename_real
    close(unit=nunit_real,iostat=ios,iomsg=iom, &
         status="keep")

    write(*,*) 'matrix2octave: closing file=',filename_imag
    close(unit=nunit_imag,iostat=ios,iomsg=iom, &
         status="keep")


500 format('# Created by matrix2octave fortran routine')
600 format('# name: matrix')
700 format('# type: matrix')
800 format('# rows: ',i10)
900 format('# columns: ',i10)
1000 format(e25.16)
    !
!!$1000 format('(',e25.16,',',e25.16')')

  end subroutine vector2octave_complex


end module matrix2octave_module
