#####-----------------------------------------

-include ../../../m_options

COMPLEX         = 0
PARALLEL        = 0
TRANSVERSE     ?= 0

SRC_PATH        = .
OBJ_ROOT        = _obj_
OBJ_PATH        = $(OBJ_ROOT)/maindir
MODULE_PATH     = $(OBJ_PATH)


# executable name
EXEC            = hyperelast_DPG_uw
EXECS           = hyperelast_DPG_uw

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)
CCINCS          = $(PROB_INCS)

###### User defined libraries
# LIBRARIES            = $(PROB_LIBS)

# Defines SRC_FILE
-include ./make.inc

SRC         = $(addprefix $(SRC_PATH)/,$(SRC_FILE))
OBJ         = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F90=.o))
# OBJ	    += $(OBJ_PATH)/main.o

ifeq ($(DEBUG),YES)
FFLAGS += -gen-interface -warn none -warn interfaces
endif
#
FC_WORK     = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

#	( cd ../.. ; make )

one : $(OBJ_PATH)/.dummy $(OBJ)

	@echo " "
	@echo "- Linking - ", $(EXEC)
	@echo "- COMPLEX, DEBUG =" $(COMPLEX), $(DEBUG)
	@echo "------------------------------ "
	$(FC_WORK) -o $(EXEC) $(OBJ) $(PROB_LIBS)

all :
	@for e in $(EXECS) ; do \
		make one EXEC=$$e ;\

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
	@rm -f $(EXECS)
	@rm -f *~

clean-pv :
	@rm -f ../outputs/paraview/*.xmf
	@rm -f ../outputs/paraview/*.h5
