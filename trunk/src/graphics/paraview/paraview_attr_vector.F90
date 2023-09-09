!-------------------------------------------------------------------------------------------
!> Purpose : driver for writing vector attribute to .h5 file
!!
!> @param[in] Id  - integer to be converted to file postfix
!> @param[in] Idx - integer identifying vector attribute
!!
!> @date Sep 2023
!-------------------------------------------------------------------------------------------
!
subroutine paraview_attr_vector(Id, Idx)
!
   use environment , only : PREFIX
   use physics     , only : NR_COMP,PHYSA
   use mpi_param   , only : RANK,ROOT
   use paraview
!
   implicit none

   integer, intent(in) :: Id
   integer, intent(in) :: Idx
!
   character(len=60) :: fname,nick
   integer           :: ic,iload,iattr,icomp,jcomp
   character(len=5)  :: postfix
   character(len=2)  :: comp,load
!
!-------------------------------------------------------------------------------------------
!
!..convert integer to string
   write(postfix,"(I5.5)") Id
!
!..decode
   iload = Idx/100
   iattr = Idx - iload*100 ; iattr=iattr/10
   icomp = Idx - iload*100 - iattr*10
   jcomp = sum(NR_COMP(1:iattr-1)) + icomp ! component index
!
!..convert integers to string
   write(comp,"(I2.2)") icomp
   write(load,"(I2.2)") iload
!
!
!  -- REAL PART --
   if (.not. PARAVIEW_COMP_REAL(jcomp)) goto 50
!
!..file name and nickname for attribute
   fname = PHYSA(iattr)//"_comp_"//comp//"_load_"//load//"_real_"
   nick  = PHYSA(iattr)//"_comp_"//comp//"_load_"//load//"_real"
!
!..write to .h5 file
   call vector2vtk("Vector",  &
      trim(PARAVIEW_DIR)//trim(PREFIX)//"_"//trim(fname)//postfix//".h5", &
      trim(nick),Idx, ic)
!
   if (RANK .ne. ROOT) goto 50
!
!..write to .xmf file (only used if XDMF/XMF format is used)
   if (.not. VIS_VTU) then
      write(PARAVIEW_IO, 1101) trim(nick)
      write(PARAVIEW_IO, 1102) ic
      write(PARAVIEW_IO, 1103) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
      write(PARAVIEW_IO, 1104)
      write(PARAVIEW_IO, 1105)
   endif
!
   50 continue
!
#if C_MODE
!
!  -- IMAGINARY PART --
   if (.not. PARAVIEW_COMP_IMAG(jcomp)) goto 70
!
!..file name and nickname for attribute
   fname = PHYSA(iattr)//"_comp_"//comp//"_load_"//load//"_imag_"
   nick  = PHYSA(iattr)//"_comp_"//comp//"_load_"//load//"_imag"
!
!..write to .h5 file (flip sign of "Idx")
   call vector2vtk("Vector",  &
      trim(PARAVIEW_DIR)//trim(PREFIX)//"_"//trim(fname)//postfix//".h5", &
      trim(nick),-abs(Idx), ic)
!
   if (RANK .ne. ROOT) goto 70
!
!..write to .xmf file (only used if XDMF/XMF format is used)
   if (.not. VIS_VTU) then
      write(PARAVIEW_IO, 1101) trim(nick)
      write(PARAVIEW_IO, 1102) ic
      write(PARAVIEW_IO, 1103) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
      write(PARAVIEW_IO, 1104)
      write(PARAVIEW_IO, 1105)
   endif
!
   70 continue
!
#endif
!
!..write to .xmf file (only used if XDMF/XMF format is used)
   if (.not. VIS_VTU) then
      1101 format("      <Attribute Name='",a,"' AttributeType='Vector' Center='Node'>")
      1102 format("        <DataItem Dimensions='",i14, " 3' NumberType='Float' Precision='4' Format='HDF'>")
      1103 format("        ",a,"vector_",a,a,".h5:",a)
      1104 format("        </DataItem>")
      1105 format("      </Attribute>")
   endif
!
   90 continue
!
end subroutine paraview_attr_vector
