subroutine vis_geom
!
      use data_structure3D , only : MAXP,NRELES
      use data_structure3D , only : NRNODS
      use environment      , only : QUIET_MODE
      use paraview         , only : PARAVIEW_DUMP_GEOM, &
                                    PARAVIEW_GEOM , &
                                    PARAVIEW_EXGEOM , &
                                    PARAVIEW_ISOGEOM
!
      implicit none
      integer :: mdle,i,ip,iq,ir
      logical :: mode_save
!
!--------------------------------------------------------------------------------
!
!     set quiet mode to TRUE
      mode_save = QUIET_MODE
      QUIET_MODE = .TRUE.
!
!     force Paraview to update geometry at each invocation
      PARAVIEW_DUMP_GEOM = .TRUE.
!
!     -- EXACT GEOMETRY --
!
!     force Paraview to use exact geometry
      PARAVIEW_GEOM = PARAVIEW_EXGEOM
!
!     dump to Paraview
      call paraview_driver
!
!
!     -- ISOPARAMETRIC GEOMETRY --
!
!     force Paraview to use isoparametric geometry
      PARAVIEW_GEOM = PARAVIEW_ISOGEOM
!
!     1st order mesh
      call set_mesh_order(1)
!
!     dump to Paraview
      call paraview_driver
!
!     subsequent orders
      do i=2,MAXP
!
!       perform uniform global p-enrichment
        call global_pref
!
!       dump to Paraview
        call paraview_driver
!
      enddo
!
!     reset quiet mode
      QUIET_MODE = mode_save
!
!
endsubroutine vis_geom
