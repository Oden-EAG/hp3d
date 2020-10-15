SRC_MODULE  = \
            FunctionalAnalysis.F90 \
            paraview_local.F90 \
            ../common/modules/common_prob_data.F90 \
            ../common/modules/geometry_transformations.F90 \
            ../common/modules/isotropic_elast_material.F90

SRC_COMMON = \
            ../common/routines/customize_geometry.F90 \
            ../common/routines/dirichlet.F90 \
            ../common/routines/paraview.F90 \
            ../common/routines/set_environment.F90 \
            ../common/routines/initialize.F90 \
            ../common/routines/finalize.F90 \
            ../common/routines/solutions.F90 \
            ../common/routines/loads.F90

SRC_ELAST  = \
            teststuff.F90 \
            elem.F90 \
            elem_residual.F90 \
            adapt_DPG.F90 \
            celem.F90 \
            exact.F90 \
            exact_error.F90 \
            exact_error_element.F90 \
            set_initial_mesh.F90 \
            soldis.F90 \
            main.F90

SRC_FILE   = $(SRC_MODULE) $(SRC_COMMON) $(SRC_ELAST)
