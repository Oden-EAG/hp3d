SRC_MODULE  = \
            ../common/modules/common_prob_data.F90 \
            ../common/modules/geometry_transformations.F90 \
            ./module/sheathed_isotropic_materials.F90 \
            ./module/primal_module.F90 \
            ./module/uweak_module.F90

SRC_COMMON = \
            ../common/routines/customize_geometry.F90 \
            ./common/dirichlet.F90 \
            ../common/routines/set_environment.F90 \
            ../common/routines/initialize.F90 \
            ../common/routines/finalize.F90 \
            ./common/solutions.F90 \
            ../common/routines/loads.F90

SRC_ELAST  = \
            ./HYBRID/elem.F90 \
            ./HYBRID/elem_residual.F90 \
            ./HYBRID/refine_DPG.F90 \
            ./HYBRID/celem.F90 \
            ./HYBRID/exact.F90 \
            ./HYBRID/exact_error.F90 \
            ./HYBRID/stress_error.F90 \
            ./HYBRID/set_initial_mesh.F90 \
            ./HYBRID/soldis.F90 \
            ./HYBRID/remove_RBM.F90 \
            ./HYBRID/main.F90 \
            ./HYBRID/elem_UW.F90 \
            ./HYBRID/elem_PRIMAL.F90 \
            ./HYBRID/paraview/paraview_custom_driver.F90 \
            ./HYBRID/paraview/paraview_attr_custom.F90 \
            ./HYBRID/paraview/soln2vtk.F90

SRC_FILE   = $(SRC_MODULE) $(SRC_COMMON) $(SRC_ELAST)
