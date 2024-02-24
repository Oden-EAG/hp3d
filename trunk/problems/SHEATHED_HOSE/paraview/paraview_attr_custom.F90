!-------------------------------------------------------------------------------------------
!> Purpose : driver for writing scalar attribute to .h5 file
!!
!> @param[in] Id   - integer to be converted to file postfix
!> @param[in Scomp - component of solution
!                       1 = sigma_rr
!                       2 = sigma_rt
!                       3 = sigma_tt
!!
!> @date Jun 16
!-------------------------------------------------------------------------------------------
!
subroutine paraview_attr_custom(Id,Scomp)
!
      use environment , only : PREFIX
      use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR
      use mpi_param   , only : RANK,ROOT
!
      implicit none
      integer, intent(in) :: Id
      integer, intent(in) :: Scomp
!
      character*8, parameter, dimension(3) :: sAttr = (/'sigma_rr','sigma_rt','sigma_tt'/)
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
      fname = trim(PREFIX)//trim(sAttr(Scomp))//'_'
      nick  = trim(PREFIX)//trim(sAttr(Scomp))
!
!     write to .h5 file
      call soln2vtk("Scalar", &
           trim(PARAVIEW_DIR)//trim(PREFIX)//"_"//trim(fname)//postfix//".h5", &
           trim(nick), "Node", Scomp, ic)
!
      if (RANK .ne. ROOT) goto 90
!
!     write to .xmf file
      write(PARAVIEW_IO, 1101) trim(nick)
      write(PARAVIEW_IO, 1102) ic
      write(PARAVIEW_IO, 1103) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
      write(PARAVIEW_IO, 1104)
      write(PARAVIEW_IO, 1105)
!
 1101 format("      <Attribute Name='",a,"' AttributeType='Scalar' Center='Node'>")
 1102 format("        <DataItem Dimensions='",i8, " 1' NumberType='Float' Precision='4' Format='HDF'>")
 1103 format("        ",a,"_",a,a,".h5:",a)
 1104 format("        </DataItem>")
 1105 format("      </Attribute>")
!
      90 continue
!
end subroutine paraview_attr_custom
