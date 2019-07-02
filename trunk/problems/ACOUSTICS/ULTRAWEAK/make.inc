SRC_MODULE  = \
            ./modules/common_prob_data.F90 \
            ./modules/geometry_transformations.F90
				

SRC_COMMON = \
            ./common/customize_geometry.F90 \
            ./common/dirichlet.F90 \
            ./common/set_environment.F90 \
            ./common/initialize.F90 \
            ./common/finalize.F90 \
            ./common/solutions.F90 \
            ./common/loads.F90

SRC_UW  = \
            elem.F90 \
            elem_fi.F90 \
            set_3Dint_fi.F90 \
            elem_residual.F90 \
            elem_residual_fi.F90 \
            refine_DPG.F90 \
            celem.F90 \
            exact.F90 \
            exact_error.F90 \
            set_initial_mesh.F90 \
            soldis.F90 \
            main.F90

SRC_SOLVER_MOD  = \
				./solver/MG/MODULE/mg_datastructure.F90 \
            ./solver/MG/MODULE/prolongation.F90 \
            ./solver/MG/MODULE/macro_grid_info.F90 \
            ./solver/MG/MODULE/constraints_info.F90 \
            ./solver/MG/MODULE/patch.F90 \
            ./solver/MG/MODULE/pcg_info.F90





SRC_SOLVER  = \
            ./solver/MG/DIRECT_SOLVE/mumps_interface.F90 \
            ./solver/MG/DIRECT_SOLVE/pardiso_interface.F90 \
            ./solver/MG/DIRECT_SOLVE/coarse_solve_mumps.F90 \
            ./solver/MG/UTILITIES/mark_masters.F90 \
            ./solver/MG/UTILITIES/find_master.F90 \
            ./solver/MG/UTILITIES/nelcon_macro.F90 \
            ./solver/MG/macro_elem_info.F90 \
            ./solver/MG/logic_macro.F90 \
            ./solver/MG/celem_sys.F90 \
            ./solver/MG/celem1.F90 \
            ./solver/MG/solout1.F90 \
            ./solver/MG/logic_nodesC.F90 \
            ./solver/MG/get_connect_infoC.F90 \
            ./solver/MG/logicC.F90 \
            ./solver/MG/coarse2macro.F90 \
            ./solver/MG/macro2coarse.F90 \
            ./solver/MG/celem_macro.F90 \
            ./solver/MG/solout_macro.F90 \
            ./solver/MG/patch_assembly.F90 \
            ./solver/MG/macro_assembly.F90 \
            ./solver/MG/debug_driver.F90 \
            ./solver/MG/DEBUG/macro_error.F90 \
            ./solver/MG/DEBUG/prolong_solution.F90 \
            ./solver/MG/DEBUG/error_contr.F90 \
            ./solver/MG/DEBUG/locate_coordC.F90 \
            ./solver/MG/DEBUG/random_refineDPG.F90 \
            ./solver/MG/mg_pcg_solve.F90 \
            ./solver/MG/mg_preconditioner.F90 \
            ./solver/MG/coarse_grid_solve.F90 \
            ./solver/MG/smoothers.F90 \
            ./solver/MG/inject_dof.F90






SRC_FILE   += $(SRC_MODULE) $(SRC_SOLVER_MOD) $(SRC_COMMON) $(SRC_UW)  $(SRC_SOLVER)
