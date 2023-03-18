# Module dependencies specified directly in makefile
#SRC_MODULE = \
#        ../COMMON_FILES/commonParam.F90 \
#        ../COMMON_FILES/laserParam.F90

SRC_COMMON = \
        ../COMMON_FILES/bessel_fn.F90 \
        ../COMMON_FILES/commonRoutines.F90 \
        ../COMMON_FILES/exact_error.F90 \
        ../COMMON_FILES/mfd_solutions.F90 \
        ../COMMON_FILES/norms.F90 \
        ../COMMON_FILES/polarization.F90 \
        ../COMMON_FILES/power.F90 \
        ../COMMON_FILES/refine_DPG.F90 \
        ../COMMON_FILES/residual.F90 \
        ../COMMON_FILES/set_environment.F90 \
        ../COMMON_FILES/thermal.F90

SRC_LASER = \
        celem.F90 \
        dirichlet.F90 \
        elem/elem.F90 \
        elem/elem_heat.F90 \
        elem/elem_maxwell.F90 \
        elem/elem_maxwell_fi_hexa.F90 \
        elem/elem_maxwell_fi_pris.F90 \
        elem/elem_maxwell_gram_hexa.F90 \
        elem/elem_maxwell_gram_pris.F90 \
        elem/elem_residual.F90 \
        elem/elem_residual_heat.F90 \
        elem/elem_residual_maxwell.F90 \
        exact.F90 \
        exec_case.F90 \
        exec_job.F90 \
        exec_job_heat.F90 \
        exec_job_nl.F90 \
        exec_job_coupled.F90 \
        getf.F90 \
        finalize.F90 \
        initialize.F90 \
        main.F90 \
        main_aux.F90 \
        my_paraview_driver.F90 \
        set_initial_mesh.F90 \
        set_3Dint_fi.F90 \
        soldis.F90
