!-------------------------------------------------------------------------------------------
!> Purpose : driver for writing geometry to .h5 file 
!!
!> @date Apr 2019
!-------------------------------------------------------------------------------------------
!
subroutine paraview_geom
!
      use environment , only : PREFIX
      use paraview    , only : PARAVIEW_IO,PARAVIEW_DUMP_GEOM,PARAVIEW_DIR
!      
      implicit none
      integer,         save :: id=-1
      character(len=5),save :: postfix
      integer,         save :: ice, icn, icp
!      
!-------------------------------------------------------------------------------------------
!
!     .h5 file is produced on 1st visit, OR as required by the user
      if (PARAVIEW_DUMP_GEOM .or. (id == -1)) then

!       increment visitation flag
        id=id+1
!      
!       convert integer to string
        write(postfix,"(I5.5)") id
!
!       produce .h5 file
        call geom2vtk("Geometry",trim(PARAVIEW_DIR)//trim(PREFIX)//"geom_"//trim(postfix)//".h5", ice, icn, icp)
!
      endif
!
!     write to .xmf file
      write(PARAVIEW_IO,1012) ice
      write(PARAVIEW_IO,1013) (ice+icn)
      write(PARAVIEW_IO,1014) trim(PREFIX), trim(postfix)
      write(PARAVIEW_IO,1015)
      write(PARAVIEW_IO,1016)
      write(PARAVIEW_IO,1017)
      write(PARAVIEW_IO,1018) icp
      write(PARAVIEW_IO,1019) trim(PREFIX), trim(postfix)
      write(PARAVIEW_IO,1020)
      write(PARAVIEW_IO,1021)
!
 1012 format("      <Topology TopologyType='Mixed' NumberOfElements='",i8,"'>")
 1013 format("        <DataItem Dimensions='",i8,"' NumberType='Int' Precision='4' Format='HDF'>")
 1014 format("        ",a,"geom_",a,".h5:/Objects")
 1015 format("        </DataItem>")
 1016 format("      </Topology>")
 1017 format("      <Geometry GeometryType='XYZ'>")
 1018 format("        <DataItem Dimensions='",i8, " 3' NumberType='Float' Precision='4' Format='HDF'>")
 1019 format("        ",a,"geom_",a,".h5:/Coords")
 1020 format("        </DataItem>")
 1021 format("      </Geometry>")
!
!
end subroutine paraview_geom
