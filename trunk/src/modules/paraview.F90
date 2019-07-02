module paraview
!
      save
! 
      character(len=128) :: FILE_VIS = '../../files/vis'
      character(len=2  ) :: VLEVEL = '2'
      character(len=128) :: PARAVIEW_DIR = 'vtk/'
      integer,parameter  :: PARAVIEW_IO = 22
      logical            :: PARAVIEW_DUMP_GEOM = .FALSE.
      logical            :: PARAVIEW_DUMP_ATTR = .FALSE.
!
!     this is matching EXGEOM flag in "control" module
      integer,parameter  :: PARAVIEW_EXGEOM  = 0
      integer,parameter  :: PARAVIEW_ISOGEOM = 1
!
!     this is reset to EXGEOM flag in "read_control" routine
      integer            :: PARAVIEW_GEOM = PARAVIEW_EXGEOM
!
!     set the domain to be output (PARAVIEW_DOMAIN = 0 means output all)
      integer            :: PARAVIEW_DOMAIN = 0
!
!
endmodule paraview
