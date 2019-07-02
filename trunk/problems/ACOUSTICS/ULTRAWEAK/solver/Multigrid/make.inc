SRC_MG_SOLVER  = \
            ./solver/Multigrid/celemMG.F90 \
            ./solver/Multigrid/celem_sysMG.F90 \
            ./solver/Multigrid/coarse_solve_mumpsMG.F90 \
            ./solver/Multigrid/elemMG.F90 \
            ./solver/Multigrid/compute_patch_mdleC.F90 \
            ./solver/Multigrid/get_sol_dofMG.F90






SRC_FILE   += $(SRC_MG_SOLVER)
