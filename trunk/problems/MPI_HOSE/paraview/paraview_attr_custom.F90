!-------------------------------------------------------------------------------------------
!> Purpose : driver for writing scalar attribute to .h5 file
!!
!> @param[in] Id   - integer to be converted to file postfix
!> @param[in Scomp - component of solution
!                       1 = sigma_rr
!                       2 = sigma_rt
!                       3 = sigma_tt
!                       4 = u_comp1
!                       5 = u_comp2
!                       6 = u_comp3
!!
!> @date Jun 16
!-------------------------------------------------------------------------------------------
!
subroutine paraview_attr_custom(Id,Scomp)
!
      use environment , only : PREFIX
      use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR
      use physics     , only : PHYSA
      use paraview    , only : PARAVIEW_DOMAIN
!
      implicit none
      integer, intent(in) :: Id
      integer, intent(in) :: Scomp
!
      character*8, parameter, dimension(4) :: sAttr = (/'sigma_rr','sigma_rt','sigma_tt','displ'/)
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
      select case(Scomp)
      case(1:3)
!
        call vis_scalar_create(    trim(PARAVIEW_DIR)//trim(PREFIX)//"scalar_"//trim(fname)//postfix//".h5", &
                               len(trim(PARAVIEW_DIR)//trim(PREFIX)//"scalar_"//trim(fname)//postfix//".h5") )
!       write to .h5 file
        call soln2vtk("Scalar", trim(fname), trim(nick), "Node", PARAVIEW_DOMAIN, Scomp, ic)
!
!       write to .xmf file
        write(PARAVIEW_IO, 1101) trim(nick)
        write(PARAVIEW_IO, 1103) ic
        write(PARAVIEW_IO, 1105) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
        write(PARAVIEW_IO, 1107)
        write(PARAVIEW_IO, 1108)
!
        call vis_scalar_release
! 
      case(4)
!
        call vis_vector_create(    trim(PARAVIEW_DIR)//trim(PREFIX)//"vector_"//trim(fname)//postfix//".h5", &
                               len(trim(PARAVIEW_DIR)//trim(PREFIX)//"vector_"//trim(fname)//postfix//".h5") )
!       write to .h5 file
        call soln2vtk("Vector", trim(fname), trim(nick), "Node", PARAVIEW_DOMAIN, Scomp, ic)
!
!       write to .xmf file
        write(PARAVIEW_IO, 1102) trim(nick)
        write(PARAVIEW_IO, 1104) ic
        write(PARAVIEW_IO, 1106) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
        write(PARAVIEW_IO, 1107)
        write(PARAVIEW_IO, 1108)    
!
        call vis_vector_release
!
      end select
!
!
 1101 format("      <Attribute Name='",a,"' AttributeType='Scalar' Center='Node'>")
 1102 format("      <Attribute Name='",a,"' AttributeType='Vector' Center='Node'>")
 1103 format("        <DataItem Dimensions='",i8, " 1' NumberType='Float' Precision='4' Format='HDF'>")
 1104 format("        <DataItem Dimensions='",i12, " 3' NumberType='Float' Precision='4' Format='HDF'>")
 1105 format("        ",a,"scalar_",a,a,".h5:",a)
 1106 format("        ",a,"vector_",a,a,".h5:",a)
 1107 format("        </DataItem>")
 1108 format("      </Attribute>")
!
!
endsubroutine paraview_attr_custom
