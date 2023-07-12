#-------------------------------------------------

-include ../../../m_options

COMPLEX   ?= 1
DEBUG     ?= 0
HP3D_USE_OPENMP    ?= 0
HP3D_USE_INTEL_MKL ?= 0

ifeq ($(COMPLEX),0)
$(error COMPLEX FLAG MUST BE 1 FOR ACOUSTICS)
endif

SRC_PATH        = .
OBJ_ROOT        = _obj_
OBJ_PATH        = $(OBJ_ROOT)/maindir
MODULE_PATH     = $(OBJ_PATH)

# Executable name
EXEC            = helm

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

# User defined libraries
LIBS            = $(PROB_LIBS)

# Defines SRC_MODULE,SRC_COMMON,SRC_LASER
-include make.inc

SRC         = $(addprefix $(SRC_PATH)/,$(SRC_FILE))
OBJ         = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F90=.o))

FC_WORK     = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

#-------------------------------------------------
.PHONY : one

one : $(OBJ_PATH)/.dummy $(OBJ)
	@echo " "
	@echo "- Linking - " $(EXEC)
	@echo "- COMPLEX = " $(COMPLEX)
	@echo "- DEBUG   = " $(DEBUG)
	@echo "------------------------------ "
	$(FC_WORK) -o $(EXEC) $(OBJ) $(LIBS) $(LDFLAGS)

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@mkdir -p $(OBJ_PATH); touch $@;
	
$(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

#-------------------------------------------------
.PHONY : print clean clean-pv

print :
	@echo $(OBJ)

clean :
	@rm -rf $(OBJ_ROOT)
	@rm -f $(EXEC)
	@rm -f *~
	@rm -f *.mod
	@rm -f *genmod.f90

clean-pv :
	@rm -f ../outputs/paraview/*
	@rm -f ../outputs/power/*
	@rm -f ../outputs/temp/*
