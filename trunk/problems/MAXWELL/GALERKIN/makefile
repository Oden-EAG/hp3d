#-------------------------------------------------

-include ../../../m_options

COMPLEX   ?= 1
DEBUG     ?= 0
HP3D_USE_OPENMP    ?= 0
HP3D_USE_INTEL_MKL ?= 0

ifeq ($(COMPLEX),0)
$(error COMPLEX FLAG MUST BE 1 FOR MAXWELL APPLICATION)
endif

SRC_PATH        = .
OBJ_ROOT        = _obj_
OBJ_PATH        = $(OBJ_ROOT)/maindir
MODULE_PATH     = $(OBJ_PATH)

# Executable name
EXEC            = maxw

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

# User defined libraries
LIBS            = $(PROB_LIBS)

# Defines SRC_MODULE,SRC_COMMON,SRC_GALERKIN
-include ./make.inc

MODULE_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_MODULE))
COMMON_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_COMMON))
GALERKIN_FILE = $(addprefix $(SRC_PATH)/,$(SRC_GALERKIN))

MODULE_OBJ    = $(addprefix $(OBJ_PATH)/, $(MODULE_FILE:.F90=.o))
COMMON_OBJ    = $(addprefix $(OBJ_PATH)/, $(COMMON_FILE:.F90=.o))
GALERKIN_OBJ  = $(addprefix $(OBJ_PATH)/, $(GALERKIN_FILE:.F90=.o))

OBJ           = $(MODULE_OBJ) $(COMMON_OBJ) $(GALERKIN_OBJ)

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

# assumes no module inter-dependencies
$(MODULE_OBJ) : $(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

# common files may depend on modules
$(COMMON_OBJ) : $(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90 $(MODULE_OBJ)
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

# Galerkin files may depend on modules
$(GALERKIN_OBJ) : $(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90 $(MODULE_OBJ)
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

#-------------------------------------------------
.PHONY : print clean clean-pv

print :
	@echo $(OBJ)

clean :
	@rm -rf $(OBJ_ROOT)
	@rm -f $(EXECS)
	@rm -f *~
	@rm -f *.mod

clean-pv :
	@rm -f ../outputs/paraview/*.xmf
	@rm -f ../outputs/paraview/*.h5
