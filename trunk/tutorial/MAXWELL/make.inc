SRC_MODULE  = \
            ./module/em.F90

SRC_CORE    = \
            celem.F90 \
            currents.F90 \
            customize_geometry.F90 \
            dirichlet.F90 \
            elem.F90 \
            initialize.F90 \
            set_initial_mesh.F90 \
            soldis.F90 

SRC_EM      = \
            elem_em.F90

SRC_FILE    = $(SRC_MODULE) $(SRC_CORE) $(SRC_EM)
