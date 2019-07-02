subroutine paraview(Fname, Title, Iex, Iap)
  use lapl
  use adaptivity, only : IPARA
  implicit none
  character(len=*), intent(in) :: Fname, Title
  integer, intent(in) :: Iex, Iap
  integer, parameter :: nfile = 18, ntmp = 19
  integer :: i
  character(len=32) :: str

  write(str,8000) (1000 + IPARA)
  open(unit=ntmp, file=trim(PREFIX)//'tmp.vtk')
  open(unit=nfile, file=trim(PREFIX)//trim(fname)//trim(str)//'.vtk')

  call mesh2vtk_geom_upscale(nfile, ISEL_PARAVIEW, 4)

  ! tmp memory flush
  call mesh2vtk_data_upscale(ntmp, Title, ISEL_PARAVIEW, &
       1.d0, mdle2vtk_data_upscale)

  if (Iex.eq.1) then
     ISEL_PARAVIEW(3) = 1
     call mesh2vtk_data_upscale(nfile, Title//'_EX', ISEL_PARAVIEW, &
          1.d0, mdle2vtk_data_upscale)

  end if

  if (Iap.eq.1) then
     ISEL_PARAVIEW(3) = 0
     call mesh2vtk_data_upscale(nfile, Title//'_AP', ISEL_PARAVIEW, &
          1.d0, mdle2vtk_data_upscale)

  end if

  close(nfile)
  close(ntmp)

8000 format(i4)
  !
end subroutine paraview

