#####-----------------------------------------

COMPLEX         = 0
PARALLEL        = 0

-include ../../../m_options

FORM            ?= PRIMAL
# FORM            = STRONG

SRC_PATH        = .
OBJ_ROOT        = _obj_
OBJ_PATH        = $(OBJ_ROOT)/maindir
MODULE_PATH     = $(OBJ_PATH)


# executable name
ifeq ($(FORM),WEAK)
EXEC            = DPG_uweak
else ifeq ($(FORM),STRONG)
EXEC            = DPG_uweak_symm
else ifeq ($(FORM),PRIMAL)
EXEC            = DPG_incomp_primal
else ifeq ($(FORM),HYBRID)
EXEC            = DPG_hybrid
else ifeq ($(FORM),PRIMALIMPROVED)
EXEC            = DPG_impr_primal
else ifeq ($(FORM),PROJ)
EXEC            = proj
endif
EXECS           = DPG_uweak DPG_uweak_symm DPG_primal DPG_impr_primal DPG_hybrid proj

# Includes and defs
CCDEFS          = $(PROB_PP_DEFS)

CCINCS          = $(PROB_INCS)


###### User defined libraries
LIBS            = $(PROB_LIBS)

# Defines SRC_FILE
ifeq ($(EXEC),DPG_uweak)
-include ./WEAK_SYMMETRY/make.inc
else ifeq ($(EXEC),DPG_uweak_symm)
-include ./STRONG_SYMMETRY/make.inc
else ifeq ($(EXEC),DPG_incomp_primal)
-include ./INCOMP_PRIMAL/make.inc
else ifeq ($(EXEC),DPG_hybrid)
-include ./HYBRID/make.inc
else ifeq ($(EXEC),DPG_impr_primal)
-include ./IMPROVED_PRIMAL/make.inc
else ifeq ($(EXEC),proj)
-include ./PROJECTION/make.inc
endif

SRC         = $(addprefix $(SRC_PATH)/,$(SRC_FILE))
OBJ         = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F90=.o))
# OBJ	    += $(OBJ_PATH)/main.o

FC_WORK     = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

#	( cd ../.. ; make )

one : $(OBJ_PATH)/.dummy $(OBJ)

	@echo " "
	@echo "- Linking - ", $(EXEC)
	@echo "- COMPLEX, PARALLEL ", $(COMPLEX) $(PARALLEL)
	@echo "------------------------------ "
	$(FC_WORK) -o $(EXEC) $(OBJ) $(LIBS)

all :
	@for e in $(EXECS) ; do \
		make one EXEC=$$e ;\
	done

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
