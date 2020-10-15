!-------------------------------------------------------------------------------------------
!> Purpose : driver for writing scalar attribute to .h5 file
!!
!> @param[in] Id   - integer to be converted to file postfix
!> @param[in Icomp - component of DISP_ATTR (see common_prob_data)
!!
!> @date Aug 2020
!-------------------------------------------------------------------------------------------
!
subroutine paraview_attr_custom(Id,Icomp)
!
  use environment , only : PREFIX
  use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR
  use physics     , only : PHYSA
  use paraview    , only : PARAVIEW_DOMAIN
  use mpi_param   , only : RANK,ROOT
  use common_prob_data
!
  implicit none
  integer, intent(in) :: Id
  integer, intent(in) :: Icomp
!
  character(len=60) :: fname,nick
  integer :: ic
  character(len=5) :: postfix
!
!-------------------------------------------------------------------------------------------
!
!     convert integer to string
  write(postfix,"(I5.5)") Id
!
!     -- REAL PART --
!
!     file name and nickname for attribute
  fname = trim(DISP_ATTR(Icomp))//"_real_"
  nick  = trim(DISP_ATTR(Icomp))//"_real"
!
!
!..write to .h5 file
  call soln2vtk("Scalar",  &
    trim(PARAVIEW_DIR)//trim(PREFIX)//"scalar_"//trim(fname)//postfix//".h5", &
    trim(nick),Icomp, ic)
!
  if (RANK .ne. ROOT) goto 50
  !
  !..write to .xmf file
  write(PARAVIEW_IO, 1101) trim(nick)
  write(PARAVIEW_IO, 1102) ic
  write(PARAVIEW_IO, 1103) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
  write(PARAVIEW_IO, 1104)
  write(PARAVIEW_IO, 1105)
  !
  50 continue
  !
  !
  if (RANK .ne. ROOT) goto 90
!
 1101 format("      <Attribute Name='",a,"' AttributeType='Scalar' Center='Node'>")
 1102 format("        <DataItem Dimensions='",i10, " 1' NumberType='Float' Precision='4' Format='HDF'>")
 1103 format("        ",a,"scalar_",a,a,".h5:",a)
 1104 format("        </DataItem>")
 1105 format("      </Attribute>")
!
   90 continue
!
!
end subroutine paraview_attr_custom
