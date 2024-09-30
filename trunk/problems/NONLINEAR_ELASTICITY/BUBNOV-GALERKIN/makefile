#####-----------------------------------------

-include ../../../m_options

COMPLEX   ?= 0
DEBUG     ?= 1
HP3D_USE_OPENMP    ?= 1
HP3D_USE_INTEL_MKL ?= 1

ifeq ($(COMPLEX),1)
$(error COMPLEX FLAG MUST BE 0 FOR HYPERELASTICITY APPLICATION)
endif

SRC_PATH        = .
OBJ_ROOT        = _obj_
OBJ_PATH        = $(OBJ_ROOT)/maindir
MODULE_PATH     = $(OBJ_PATH)

# executable name
EXEC            = hyperelast_BG

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

# User defined libraries
LIBS            = $(PROB_LIBS)

# Defines SRC_MODULE,SRC_COMMON,SRC_HYPER
-include ./make.inc

MODULE_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_MODULE))
COMMON_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_COMMON))
HYPER_FILE    = $(addprefix $(SRC_PATH)/,$(SRC_HYPER))

MODULE_OBJ    = $(addprefix $(OBJ_PATH)/, $(MODULE_FILE:.F90=.o))
COMMON_OBJ    = $(addprefix $(OBJ_PATH)/, $(COMMON_FILE:.F90=.o))
HYPER_OBJ     = $(addprefix $(OBJ_PATH)/, $(HYPER_FILE:.F90=.o))

OBJ           = $(MODULE_OBJ) $(COMMON_OBJ) $(HYPER_OBJ)

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

# elast files may depend on modules
$(HYPER_OBJ) : $(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90 $(MODULE_OBJ)
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

clean-pv :
	@rm -f output/figures/*.xmf
	@rm -f output/figures/*.h5
