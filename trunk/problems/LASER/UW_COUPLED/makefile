#####-----------------------------------------

-include ../../../m_options

COMPLEX   ?= 1
DEBUG     ?= 0

ifeq ($(COMPLEX),0)
$(error COMPLEX FLAG MUST BE 1 FOR LASER APPLICATION)
endif

ifeq ($(COMPILER),INTEL)
ifeq ($(DEBUG),1)
FFLAGS += -gen-interface -warn none -warn interfaces
endif
endif

SRC_PATH        = .
#OBJ_PATH        = _obj_
OBJ_PATH        = COMMON_FILES
MODULE_PATH     = $(OBJ_PATH)

# Executable name
EXEC            = uwLaser

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

# User defined libraries
LIBS            = $(PROB_LIBS)

include make.inc

SRC         = $(addprefix $(SRC_PATH)/,$(SRC_FILE))

SRC_FILE_F  = $(filter %.F, $(SRC_FILE))
SRC_FILE_F90= $(filter %.F90, $(SRC_FILE))

OBJ_F90     = $(addprefix $(OBJ_PATH)/, $(SRC_FILE_F90:.F90=.F90.o))
OBJ_F       = $(addprefix $(OBJ_PATH)/, $(SRC_FILE_F:.F=.F.o))

OBJ         = $(OBJ_F) $(OBJ_F90)


FC_WORK     = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

one : $(OBJ_PATH)/.dummy $(OBJ)
	@echo " "
	@echo "- Linking , " $(EXEC)
	@echo "- COMPLEX = " $(COMPLEX)
	@echo "- DEBUG   = " $(DEBUG)
	@echo "------------------------------ "
	$(FC_WORK) -o $(EXEC) $(OBJ) $(LIBS) $(LDFLAGS)

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_PATH); touch $@; \
	fi


$(OBJ_PATH)/%.F.o : $(SRC_PATH)/%.F
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<


$(OBJ_PATH)/%.F90.o : $(SRC_PATH)/%.F90
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
	@rm -f ../outputs/paraview/*
	@rm -f ../outputs/power/*
	@rm -f ../outputs/temp/*

