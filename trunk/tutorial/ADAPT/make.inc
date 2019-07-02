SRC_MODULE  = \
            ./module/environment.F90 \
            ./module/lapl.F90 \
            ./module/adaptivity.F90

SRC_CORE    = \
            celem.F90 \
            customize_geometry.F90 \
            dirichlet.F90 \
            elem.F90 \
            exact.F90 \
            exact_error.F90 \
            exact_error_element.F90 \
            initialize.F90 \
            loads.F90 \
            paraview.F90 \
            set_initial_mesh.F90 \
            soldis.F90 


SRC_LAPL   = \
            elem_lapl.F90

SRC_FILE    = $(SRC_MODULE) $(SRC_CORE) $(SRC_LAPL)
