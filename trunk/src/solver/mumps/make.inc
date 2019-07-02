OBJMUMPS_INTERF = $(SRC_DIR)/mumps/mumps_data.F \
                  $(SRC_DIR)/mumps/initialize_mumps.F

ifeq ($(PARALLEL),1)
    OBJMUMPS_INTERF += $(SRC_DIR)/mumps/mumps_solve_ASSEMBLED_PARALLEL.F 
else
    OBJMUMPS_INTERF += $(SRC_DIR)/mumps/mumps_solve_ASSEMBLED_SEQUENTIAL.F
endif

objmumps = $(OBJMUMPS_INTERF)

