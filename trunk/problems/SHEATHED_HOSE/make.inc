SRC_MODULE  = \
            ./module/common_prob_data.F90 \
            ./module/sheathed_isotropic_materials.F90 \
            ./module/primal_module.F90 \
            ./module/uweak_module.F90 \
            ./paraview/paraview_attr_custom.F90 \
            ./paraview/paraview_custom_driver.F90 \
            ./paraview/soln2vtk.F90

SRC_COMMON = \
            ./common/dirichlet.F90 \
            ./common/set_environment.F90 \
            ./common/initialize.F90 \
            ./common/finalize.F90 \
            ./common/solution.F90 \
            ./common/residual.F90 \
            ./common/uniform_href.F90 \
            ./common/getf.F90

SRC_ELAST  = \
            elem.F90 \
            elem_residual.F90 \
            exec_case.F90 \
            exec_job.F90 \
            celem.F90 \
            exact.F90 \
            exact_error.F90 \
            stress_error.F90 \
            set_initial_mesh.F90 \
            soldis.F90 \
            remove_RBM.F90 \
            main.F90 \
            elem_UW.F90 \
            elem_PRIMAL.F90 \
            my_paraview_driver.F90
