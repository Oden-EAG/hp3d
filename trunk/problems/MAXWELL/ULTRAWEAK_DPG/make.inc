SRC_MODULE  = \
        modules/commonParam.F90

SRC_COMMON  = \
        common/commonRoutines.F90 \
        common/exact_error.F90 \
        common/mfd_solutions.F90 \
        common/residual.F90 \
        common/set_environment.F90

SRC_UW      = \
        celem.F90 \
        dirichlet.F90 \
        elem/elem.F90 \
        elem/elem_maxwell.F90 \
        elem/elem_residual.F90 \
        elem/elem_residual_maxwell.F90 \
        exact.F90 \
        exec_case.F90 \
        exec_job.F90 \
        getf.F90 \
        finalize.F90 \
        initialize.F90 \
        main.F90 \
        my_paraview_driver.F90 \
        refine_DPG.F90 \
        set_initial_mesh.F90 \
        soldis.F90
