#####-----------------------------------------

-include ../../../m_options

COMPLEX        ?= 0
DEBUG_MODE     ?= 0
HP3D_USE_OPENMP    ?= 0
HP3D_USE_INTEL_MKL ?= 0

ifeq ($(COMPLEX),1)
$(error COMPLEX FLAG MUST BE 0 FOR POISSON APPLICATION)
endif

ifeq ($(COMPILER),INTEL)
ifeq ($(DEBUG),1)
FFLAGS += -gen-interface -warn none -warn interfaces
endif
endif

SRC_PATH        = .
OBJ_ROOT        = _obj_
OBJ_PATH        = $(OBJ_ROOT)/maindir
MODULE_PATH     = $(OBJ_PATH)

# Executable name
EXEC            = pois

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

###### User defined libraries
LIBS            = $(PROB_LIBS)

# Defines SRC_FILE
-include ./make.inc

SRC         = $(addprefix $(SRC_PATH)/,$(SRC_FILE))
OBJ         = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F90=.o))
# OBJ	    += $(OBJ_PATH)/main.o

FC_WORK     = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

one : $(OBJ_PATH)/.dummy $(OBJ)
	@echo " "
	@echo "- Linking - ", $(EXEC)
	@echo "- COMPLEX = " $(COMPLEX)
	@echo "- DEBUG   = " $(DEBUG)
	@echo "------------------------------ "
	$(FC_WORK) -o $(EXEC) $(OBJ) $(LIBS) $(LDFLAGS)

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_ROOT); mkdir $(OBJ_PATH); touch $@; \
	fi

$(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

print :
	@echo $(OBJ)

clean :
	@rm -rf $(OBJ_ROOT)
	@rm -f $(EXEC)
	@rm -f *~
	@rm -f *.mod

clean-pv :
	@rm -f ../outputs/paraview/*.xmf
	@rm -f ../outputs/paraview/*.h5
