#-------------------------------------------------

-include ../../m_options

COMPLEX   ?= 0
DEBUG     ?= 0
HP3D_USE_OPENMP    ?= 0
HP3D_USE_INTEL_MKL ?= 0

ifeq ($(COMPLEX),1)
$(error COMPLEX FLAG MUST BE 0 FOR LINEAR ELASTICITY APPLICATION)
endif

SRC_PATH        = .
OBJ_ROOT        = _obj_
OBJ_PATH        = $(OBJ_ROOT)/maindir
MODULE_PATH     = $(OBJ_PATH)

# Executable name
EXEC            = hose

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

# User defined libraries
LIBS            = $(PROB_LIBS)

# Defines SRC_MODULE,SRC_COMMON,SRC_ELAST
-include ./make.inc

MODULE_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_MODULE))
COMMON_FILE   = $(addprefix $(SRC_PATH)/,$(SRC_COMMON))
ELAST_FILE    = $(addprefix $(SRC_PATH)/,$(SRC_ELAST))

MODULE_OBJ    = $(addprefix $(OBJ_PATH)/, $(MODULE_FILE:.F90=.o))
COMMON_OBJ    = $(addprefix $(OBJ_PATH)/, $(COMMON_FILE:.F90=.o))
ELAST_OBJ     = $(addprefix $(OBJ_PATH)/, $(ELAST_FILE:.F90=.o))

OBJ           = $(MODULE_OBJ) $(COMMON_OBJ) $(ELAST_OBJ)

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
$(ELAST_OBJ) : $(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90 $(MODULE_OBJ)
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
	@rm -f ../outputs/paraview/*.vtu
	@rm -f ../outputs/paraview/*.pvd
