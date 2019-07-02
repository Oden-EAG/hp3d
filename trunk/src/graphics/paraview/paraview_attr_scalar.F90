!-------------------------------------------------------------------------------------------
!> Purpose : driver for writing scalar attribute to .h5 file
!!
!> @param[in] Id  - integer to be converted to file postfix
!> @param[in] Idx - index identifying scalar attribute
!!
!> @date Apr 2019
!-------------------------------------------------------------------------------------------
!
subroutine paraview_attr_scalar(Id, Idx)
!
      use environment , only : PREFIX
      use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR
      use physics     , only : PHYSA 
      use paraview    , only : PARAVIEW_DOMAIN
!
      implicit none
      integer, intent(in) :: Id
      integer, intent(in) :: Idx
!
      character(len=60) :: fname,nick
      integer :: ic,iload,iattr,icomp
      character(len=5) :: postfix
      character(len=2) :: comp,load
!
!-------------------------------------------------------------------------------------------
!
!     convert integer to string
      write(postfix,"(I5.5)") Id
!
!     decode
      iload = Idx/100
      iattr = Idx - iload*100 ; iattr=iattr/10
      icomp = Idx - iload*100 - iattr*10
!
!     convert integers to string
      write(comp,"(I2.2)") icomp
      write(load,"(I2.2)") iload
!
!
!     -- REAL PART --
!
!     file name and nickname for attribute 
      fname = PHYSA(iattr)//"_comp_"//comp//"_load_"//load//"_real_"
      nick  = PHYSA(iattr)//"_comp_"//comp//"_load_"//load//"_real"
!
      call vis_scalar_create(    trim(PARAVIEW_DIR)//trim(PREFIX)//"scalar_"//trim(fname)//postfix//".h5", &
                             len(trim(PARAVIEW_DIR)//trim(PREFIX)//"scalar_"//trim(fname)//postfix//".h5") )
!
!     write to .h5 file
      call scalar2vtk("Scalar", trim(fname), trim(nick), "Node", PARAVIEW_DOMAIN, Idx, ic); 
!
!     write to .xmf file
      write(PARAVIEW_IO, 1101) trim(nick)
      write(PARAVIEW_IO, 1102) ic
      write(PARAVIEW_IO, 1103) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
      write(PARAVIEW_IO, 1104)
      write(PARAVIEW_IO, 1105)
!
      call vis_scalar_release
!
!
#if C_MODE
!
!     -- IMAGINARY PART --
!
!     file name and nickname for attribute 
      fname = PHYSA(iattr)//"_comp_"//comp//"_load_"//load//"_imag_"
      nick  = PHYSA(iattr)//"_comp_"//comp//"_load_"//load//"_imag"
!
      call vis_scalar_create(    trim(PARAVIEW_DIR)//trim(PREFIX)//"scalar_"//trim(fname)//postfix//".h5", &
                             len(trim(PARAVIEW_DIR)//trim(PREFIX)//"scalar_"//trim(fname)//postfix//".h5") )
!
!     write to .h5 file (flip sign of "Idx")
      call scalar2vtk("Scalar", trim(fname), trim(nick), "Node", PARAVIEW_DOMAIN, -iabs(Idx), ic); 
!
!     write to .xmf file
      write(PARAVIEW_IO, 1101) trim(nick)
      write(PARAVIEW_IO, 1102) ic
      write(PARAVIEW_IO, 1103) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
      write(PARAVIEW_IO, 1104)
      write(PARAVIEW_IO, 1105)
!
      call vis_scalar_release
!
#endif
   90 continue
!
 1101 format("      <Attribute Name='",a,"' AttributeType='Scalar' Center='Node'>")
 1102 format("        <DataItem Dimensions='",i8, " 1' NumberType='Float' Precision='4' Format='HDF'>")
 1103 format("        ",a,"scalar_",a,a,".h5:",a)
 1104 format("        </DataItem>")
 1105 format("      </Attribute>")
!
!
end subroutine paraview_attr_scalar
