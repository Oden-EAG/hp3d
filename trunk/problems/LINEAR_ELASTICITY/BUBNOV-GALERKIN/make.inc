SRC_MODULE  = \
            ../common/modules/common_prob_data.F90 \
            ../common/modules/geometry_transformations.F90 \
            ../common/modules/isotropic_elast_material.F90 \
            ../common/modules/transverse_isotropic_elast_material.F90

SRC_COMMON = \
            ../common/routines/uniform_href.F90 \
            ../common/routines/customize_geometry.F90 \
            ../common/routines/dirichlet.F90 \
            ../common/routines/set_environment_galerkin.F90 \
            ../common/routines/initialize.F90 \
            ../common/routines/finalize.F90 \
            ../common/routines/solutions.F90 \
            ../common/routines/my_paraview_driver.F90 \
            ../common/routines/loads.F90

SRC_ELAST  = \
            elem.F90 \
            celem.F90 \
            exact.F90 \
            exact_error.F90 \
            set_initial_mesh.F90 \
            exec_case.F90 \
            exec_job.F90 \
            soldis.F90 \
            main.F90

SRC_FILE   = $(SRC_MODULE) $(SRC_COMMON) $(SRC_ELAST)
