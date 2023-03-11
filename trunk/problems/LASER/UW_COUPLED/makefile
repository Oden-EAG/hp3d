#-------------------------------------------------

-include ../../../m_options

COMPLEX   ?= 1
DEBUG     ?= 0
HP3D_USE_OPENMP    ?= 0
HP3D_USE_INTEL_MKL ?= 0

ifeq ($(COMPLEX),0)
$(error COMPLEX FLAG MUST BE 1 FOR LASER APPLICATION)
endif

SRC_PATH        = .
OBJ_ROOT        = _obj_
OBJ_PATH        = $(OBJ_ROOT)/maindir
MODULE_PATH     = $(OBJ_PATH)

# Executable name
EXEC            = uwLaser

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

# User defined libraries
LIBS            = $(PROB_LIBS)

# Defines SRC_MODULE,SRC_COMMON,SRC_LASER
-include make.inc

# Modules (resolve dependencies between modules)
commonParam_SRC = $(SRC_PATH)/../COMMON_FILES/commonParam.F90
laserParam_SRC  = $(SRC_PATH)/../COMMON_FILES/laserParam.F90

commonParam_OBJ = $(OBJ_PATH)/$(SRC_PATH)/../COMMON_FILES/commonParam.o
laserParam_OBJ  = $(OBJ_PATH)/$(SRC_PATH)/../COMMON_FILES/laserParam.o

MODULE_OBJ    = $(commonParam_OBJ) $(laserParam_OBJ)

COMMON_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_COMMON))
LASER_FILE    = $(addprefix $(SRC_PATH)/,$(SRC_LASER))

COMMON_OBJ    = $(addprefix $(OBJ_PATH)/, $(COMMON_FILE:.F90=.o))
LASER_OBJ     = $(addprefix $(OBJ_PATH)/, $(LASER_FILE:.F90=.o))

OBJ           = $(MODULE_OBJ) $(COMMON_OBJ) $(LASER_OBJ)

FC_WORK       = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

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

### Module dependencies ###
$(commonParam_OBJ) : $(commonParam_SRC)
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

$(laserParam_OBJ) : $(laserParam_SRC) $(commonParam_OBJ)
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<
### End Module dependencies ###

# common files may depend on modules
$(COMMON_OBJ) : $(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90 $(MODULE_OBJ)
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

# laser files may depend on modules
$(LASER_OBJ) : $(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90 $(MODULE_OBJ)
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
