SRC_MODULE  = \
            module/proj.F90

SRC_CORE    = \
            celem.F90 \
            customize_geometry.F90 \
            dirichlet.F90 \
            elem.F90 \
            exact.F90 \
            exact_error.F90 \
            exact_error_element.F90 \
            geometry_error.F90 \
            geometry_error_element.F90 \
            initialize.F90 \
            set_initial_mesh.F90 \
            soldis.F90 

SRC_PROJ   = \
            elem_proj.F90 

SRC_FILE    = $(SRC_MODULE) $(SRC_CORE) $(SRC_PROJ)
