SRC_MODULE  = \
            ./module/dome.F90

SRC_CORE    = \
            celem.F90 \
	    customize_geometry.F90 \
            dirichlet.F90 \
            elem.F90 \
            initialize.F90 \
            loads.F90 \
            set_initial_mesh.F90 \
            soldis.F90 

SRC_ELAST   = \
            elem_elast.F90

SRC_FILE    = $(SRC_MODULE) $(SRC_CORE) $(SRC_ELAST)
