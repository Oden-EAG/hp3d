!-------------------------------------------------------------------------------------------
!> Purpose : driver for writing geometry to .h5 file
!!
!> @date Mar 2023
!-------------------------------------------------------------------------------------------
!
subroutine paraview_geom
!
   use environment     , only : PREFIX
   use paraview        , only : PARAVIEW_IO, PARAVIEW_DUMP_GEOM,   &
                                PARAVIEW_DIR, SECOND_ORDER_VIS
   use mpi_param       , only : RANK, ROOT
   use data_structure3D, only : NODES, ELEM_ORDER
   use node_types
!
   implicit none
   integer,         save :: id=-1
   character(len=5),save :: postfix
   integer,         save :: ntype, ice, icn, icp, ico, mdle
!
!-------------------------------------------------------------------------------------------
!
!...h5 file is produced on 1st visit, OR as required by the user
   if (PARAVIEW_DUMP_GEOM .or. (id == -1)) then
!
!  ...increment visitation flag
      id=id+1
!
!  ...convert integer to string
      write(postfix,"(I5.5)") id
!
!  ...produce .h5 file
      call geom2vtk("Geometry",trim(PARAVIEW_DIR)//trim(PREFIX)//"geom_"//trim(postfix)//".h5", ice,icn,icp)
!
   endif
!
   if (RANK .ne. ROOT) goto 90
!
   if (SECOND_ORDER_VIS) then
      ico = icn
!
      mdle = ELEM_ORDER(1)
      ntype = NODES(mdle)%ntype
!
      select case(ntype)
      case(MDLB)
         write(PARAVIEW_IO,1011) "'HEXAHEDRON_27'", ice, 27
      case(MDLP)
         write(PARAVIEW_IO,1011) "'WEDGE_18'", ice, 18
      case(MDLN)
         write(PARAVIEW_IO,1011) "'TETRAHEDRON_10'", ice, 10
      case default
         write(*,*) 'paraview_geom: unrecognized element type: ', S_Type(ntype)
         stop 1
      end select
   else
      ico = (ice+icn)
      write(PARAVIEW_IO,1012) "'Mixed'", ice
   endif
!
!..write to .xmf file
   write(PARAVIEW_IO,1013) ico
   write(PARAVIEW_IO,1014) trim(PREFIX), trim(postfix)
   write(PARAVIEW_IO,1015)
   write(PARAVIEW_IO,1016)
   write(PARAVIEW_IO,1017)
   write(PARAVIEW_IO,1018) icp
   write(PARAVIEW_IO,1019) trim(PREFIX), trim(postfix)
   write(PARAVIEW_IO,1020)
   write(PARAVIEW_IO,1021)
!
 1011 format("      <Topology TopologyType=",a," NumberOfElements='",i8,"' NodesPerElement='",i8,"'>")
 1012 format("      <Topology TopologyType=",a," NumberOfElements='",i8,"'>")
 1013 format("        <DataItem Dimensions='",i12,"' NumberType='Int' Precision='4' Format='HDF'>")
 1014 format("        ",a,"geom_",a,".h5:/Objects")
 1015 format("        </DataItem>")
 1016 format("      </Topology>")
 1017 format("      <Geometry GeometryType='XYZ'>")
 1018 format("        <DataItem Dimensions='",i10, " 3' NumberType='Float' Precision='4' Format='HDF'>")
 1019 format("        ",a,"geom_",a,".h5:/Coords")
 1020 format("        </DataItem>")
 1021 format("      </Geometry>")
!
   90 continue
!
end subroutine paraview_geom
