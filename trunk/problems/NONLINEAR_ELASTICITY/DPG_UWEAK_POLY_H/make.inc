SRC_COMMON = \
            ../common/routines/set_environment_poly.F90 \
            ../common/routines/initialize_polydpg.F90 \
            ./find_material.F90 \
            ../common/routines/solutions.F90 \
            ../common/routines/loads.F90 \
            exact_error.F90 \
            generate_visual_poly_ct.F90 \
            Dirichlet_f_proj.F90 \
            soldis_poly.F90 \
            ./verify/poly_geom_util_verify.F90 \
            ./verify/initialize_verify.F90 \
            ./verify/finalize_verify.F90 \
            set_bc_surf.F90

SRC_MODULE = \
            ../common/modules/common_prob_data.F90 \
            ../common/modules/isotropic_elast_material_compo.F90 \
            metis_tests.F90 \
            uweak_module_poly.F90

SRC_ELAST  = \
            elem_poly_hni.F90 \
            elem_residual_hni.F90 \
            exact_poly.F90 \
            correct_dp.F90 \
            main.F90

SRC_FILE   = $(SRC_MODULE) $(SRC_COMMON) $(SRC_ELAST)
