SRC_MODULE  = \
            ./module/lapl.F90

SRC_CORE    = \
            celem.F90 \
            compute_error.F90 \
            customize_geometry.F90 \
            dirichlet.F90 \
            elem.F90 \
            error_h1.F90 \
            initialize.F90 \
            loads.F90 \
            set_initial_mesh.F90 \
            soldis.F90 \
            write_vtk_data.F90

SRC_LAPL   = \
            elem_lapl.F90

SRC_FILE    = $(SRC_MODULE) $(SRC_CORE) $(SRC_LAPL)
