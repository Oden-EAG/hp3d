SRC_MODULE = \
            ./modules/common_prob_data.F90

SRC_COMMON = \
            ./common/dirichlet.F90 \
            ./common/exact.F90 \
            ./common/exact_error.F90 \
            ./common/finalize.F90 \
            ./common/getf.F90 \
            ./common/initialize.F90 \
            ./common/residual.F90 \
            ./common/solution.F90 \
            ./common/set_environment.F90 \
            ./common/uniform_href.F90 \
            ./common/ficheracornerdirichlet.F90 \
            ./common/refine_DPG.F90
            

SRC_UW = \
            celem.F90 \
            exec_case.F90 \
            exec_job.F90 \
            exec_job_adap_ref.F90 \
            elem.F90 \
            elem_residual.F90 \
            main.F90 \
            my_paraview_driver.F90 \
            set_initial_mesh.F90 \
            soldis.F90