SRC_MODULE  = \
            modules/common_prob_data.F90 \
            modules/isotropic_elast_material.F90 \

SRC_COMMON = \
            common/set_environment.F90 \
            common/initialize.F90 \
            common/finalize.F90 \
            common/solutions.F90 \
            common/loads.F90

SRC_ELAST  = \
            celem.F90 \
            dirichlet.F90 \
            elem.F90 \
            elem_residual.F90 \
            exact.F90 \
            exact_error.F90 \
            exec_case.F90 \
            exec_job.F90 \
            main.F90 \
            my_paraview_driver.F90 \
            refine_DPG.F90 \
            set_initial_mesh.F90 \
            soldis.F90
