SRC_MODULE  = \
            ./modules/common_prob_data.F90 \
            ./modules/geometry_transformations.F90 \

SRC_COMMON = \
            ./common/customize_geometry.F90 \
            ./common/dirichlet.F90 \
            ./common/set_environment.F90 \
            ./common/initialize.F90 \
            ./common/finalize.F90 \
            ./common/solutions.F90 \
            ./common/loads.F90

SRC_PROB    = \
            elem.F90 \
            elem_residual.F90 \
            refine_DPG.F90 \
            celem.F90 \
            exact.F90 \
            exact_error.F90 \
            set_initial_mesh.F90 \
            soldis.F90 \
            main.F90

SRC_FILE   = $(SRC_MODULE) $(SRC_COMMON) $(SRC_PROB)
