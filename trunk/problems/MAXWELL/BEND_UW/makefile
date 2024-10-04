#####-----------------------------------------

-include ../../../m_options

HP3D_COMPLEX   ?= 1
HP3D_DEBUG     ?= 0

ifeq ($(HP3D_COMPLEX),0)
$(error HP3D_COMPLEX FLAG MUST BE 1 FOR MAXWELL)
endif

ifeq ($(COMPILER),INTEL)
ifeq ($(HP3D_DEBUG),1)
FFLAGS += -gen-interface -warn none -warn interfaces
endif
endif

SRC_PATH        = .
OBJ_PATH        = _obj_
MODULE_PATH     = $(OBJ_PATH)

# Executable name
EXEC            = bending

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

# User defined libraries
LIBS            = $(PROB_LIBS)

include make.inc

MODULE_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_MODULE))
COMMON_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_COMMON))
UW_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_UW))

MODULE_OBJ    = $(addprefix $(OBJ_PATH)/, $(MODULE_FILE:.F90=.o))
COMMON_OBJ    = $(addprefix $(OBJ_PATH)/, $(COMMON_FILE:.F90=.o))
UW_OBJ    = $(addprefix $(OBJ_PATH)/, $(UW_FILE:.F90=.o))

OBJ           = $(MODULE_OBJ) $(COMMON_OBJ) $(UW_OBJ)

FC_WORK       = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

one : $(OBJ_PATH)/.dummy $(OBJ)
	@echo " "
	@echo "- Linking , " $(EXEC)
	@echo "- HP3D_COMPLEX = " $(HP3D_COMPLEX)
	@echo "- HP3D_DEBUG   = " $(HP3D_DEBUG)
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

# UW files may depend on modules
$(UW_OBJ) : $(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90 $(MODULE_OBJ)
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<


clean :
	@rm -rf $(OBJ_PATH)
	@rm -f $(EXEC)
	@rm -f *~
	@rm -f *.mod
	@rm -f *genmod.f90

clean-pv :
	@rm -f ../outputs/paraview/*.h5
	@rm -f ../outputs/paraview/*.xmf
	@rm -f ../outputs/paraview/*.vtu
	@rm -f ../outputs/paraview/*.pvd
