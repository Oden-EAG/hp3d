module PML
  !  ...PML PARAMETERS
  real(8) :: RPML_MIN, RPML_MAX
  !  ...PML domain, sphere
  integer :: PML_DOMAIN, PML_SPHERE
contains

  subroutine disp_pml_data
    write(*, 6000)
    write(*, 6001) PML_DOMAIN, PML_SPHERE
    write(*, 6002) RPML_MIN, RPML_MAX
6000 format('-------------------------')
6001 format(' PML_DOMAIN, SPHERE  :   ', 2(i3,2x))
6002 format(' PML_MIN, MAX        :   ', 2(e12.5,2x))
  end subroutine disp_pml_data

  !  ...dump out the PML parameters
  subroutine dumpout_PML
    integer, parameter :: ndump=31
    open(unit=ndump,file='files/dumpPML',form='formatted',access='sequential',status='unknown')
    write(ndump,*) RPML_MIN,RPML_MAX
    close(ndump)
  end subroutine dumpout_PML
  subroutine dumpin_PML
    integer, parameter :: ndump=31
    open(unit=ndump,file='files/dumpGMP', form='formatted',access='sequential',status='unknown')
    read(ndump,*)  RPML_MIN,RPML_MAX
    close(ndump)
  end subroutine dumpin_PML
end module PML

