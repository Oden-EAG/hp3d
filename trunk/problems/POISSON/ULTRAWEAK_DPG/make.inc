SRC_MODULE = \
            ./modules/common_prob_data.F90

SRC_COMMON = \
            ./common/customize_geometry.F90 \
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
            ./common/refine_DPG.F90 \
            ./common/HpAdapt.F90 \
            ./common/project_p.F90  \
            ./common/extraction_vector.F90  \
            ./common/pbisolver3.F90 \
            ./common/Finegrid_padap.F90 \
            ./common/fine_to_coarse_gp_map.F90  \
            ./common/fine_to_coarse_projection_error.F90 \
            ./common/project_h.F90  \
            ./common/elem_proj_h.F90 \
            ./common/mesh_copy.F90  \
            ./common/poly_adap_subson.F90   \
            ./common/coeff_polynomial_aniso.F90

SRC_PRIMAL = \
            celem.F90 \
            exec_case.F90 \
            exec_job.F90 \
            elem.F90 \
            elem_residual.F90 \
            main.F90 \
            my_paraview_driver.F90 \
            set_initial_mesh.F90 \
            soldis.F90

SRC_FILE   += $(SRC_MODULE) $(SRC_COMMON) $(SRC_PRIMAL)
