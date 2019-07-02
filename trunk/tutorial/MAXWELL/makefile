#####-----------------------------------------

COMPLEX         = 1
PARALLEL        = 0


-include ../../m_options

SRC_PATH        = .
OBJ_PATH        = _obj_
MODULE_PATH     = $(OBJ_PATH)


# executable name
EXEC            = em

# Includes and defs
CCDEFS          = -D C_MODE=$(COMPLEX) \
                  -D PARALLEL_MODE=$(PARALLEL)

CCINCS          = -I$(HP3D_PATH)/common/hp3d \
                  -I$(HP3D_PATH)/common/graphics \
                  -I$(HP3D_PATH)/module \
		  -I$(MUMPS_MPI_INC) \
		  -I./$(MODULE_PATH)

###### User defined libraries
LIBS            = $(HP3D_LIB) $(HP3D_COMMON) $(HP3D_GMP) $(HP3D_COMMON) \
		  $(HP3D_LIB) \
                  $(UHM_LIB) $(FLA_LIB) \
		  $(MUMPS_LIB) $(MUMPS_MPI_LIB) \
		  $(WSMP_LIB) $(PARDISO_LIB) \
		  $(METIS_LIB) \
		  $(MKL_LIB) $(LAPACK_LIB) $(BLAS_LIB) \
		  $(X_LIB) $(PTHREAD_LIB)

ifdef CUDA_PATH
LIBS           += $(CUDA_LIB)
endif


-include make.inc

SRC         = $(addprefix $(SRC_PATH)/,$(SRC_FILE))
OBJ         = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F90=.o)) 
OBJ	    += $(OBJ_PATH)/main_$(EXEC).o

FC_WORK     = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

#	( cd ../.. ; make )

one : $(OBJ_PATH)/.dummy $(OBJ)

	@echo " "
	@echo "- Linking - ", $(EXEC)
	@echo "- COMPLEX, PARALLEL ", $(COMPLEX) $(PARALLEL)
	@echo "------------------------------ "
	$(FC_WORK) -o $(EXEC) $(OBJ) $(LIBS)  $(LDFLAGS)

all :
	@for e in $(EXECS) ; do \
		make one EXEC=$$e ;\
	done

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH) 
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_PATH); touch $@; \
	fi

$(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

print : 
	@echo $(OBJ)

clean :
	@rm -rf $(OBJ_PATH)
	@rm -f $(EXECS)
	@rm -f *~
