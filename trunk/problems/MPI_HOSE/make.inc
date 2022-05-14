SRC_MODULE  = \
            ./module/common_prob_data.F90 \
            ./module/geometry_transformations.F90 \
            ./module/sheathed_isotropic_materials.F90 \
            ./module/primal_module.F90 \
            ./module/uweak_module.F90

SRC_COMMON = \
            ./common/routines/customize_geometry.F90 \
            ./common/dirichlet.F90 \
            ./common/routines/set_environment.F90 \
            ./common/routines/initialize.F90 \
            ./common/routines/finalize.F90 \
            ./common/solutions.F90 \
            ./common/routines/getf.F90

SRC_ELAST  = \
            elem.F90 \
            elem_residual.F90 \
            refine_DPG.F90 \
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

SRC_FILE   = $(SRC_MODULE) $(SRC_COMMON) $(SRC_ELAST)
