#-------------------------------------------------
# - Library path -

-include ./m_options

# Set variables if not set by m_options
COMPLEX    ?= 0
DEBUG      ?= 0

HP3D_USE_OPENMP    ?= 0
HP3D_USE_INTEL_MKL ?= 0

####################################################

HP3D_LIBNAME   = libhp3d.a
FRONTAL_NAME   = libfrontal.a
HP3D_COMNAME   = libhp3d_common.a
HP3D_GMPNAME   = libhp3d_gmp.a

#-------------------------------------------------
# Complex version of library
ifeq ($(COMPLEX),1)
INSTALL_PATH   = ./complex
CCINCS         = -I./$(SRC_PATH)/common/frontal/complex
FRONTAL_PATH   = ./$(SRC_PATH)/solver/frontal/base/complex
endif
#-------------------------------------------------
# Real version of library
ifeq ($(COMPLEX),0)
INSTALL_PATH   = ./real
CCINCS         = -I./$(SRC_PATH)/common/frontal/real
FRONTAL_PATH   = ./$(SRC_PATH)/solver/frontal/base/real
endif
#-------------------------------------------------
# - Compiler -
# -- Preprocessor options -- 
# Files with ending in .F and .F90 (not .f and .f90)
# are preprocessed according to certain directives
CCDEFS    = -D"C_MODE=$(COMPLEX)" \
            -D"DEBUG_MODE=$(DEBUG)" \
            -D"HP3D_USE_OPENMP=$(HP3D_USE_OPENMP)" \
            -D"HP3D_USE_INTEL_MKL=$(HP3D_USE_INTEL_MKL)"

# -- Includes --
CCINCS   += -I./$(SRC_PATH) \
            -I./$(SRC_PATH)/common/graphics \
            -I./$(SRC_PATH)/common/hp3d \
            -I./$(MODULE_PATH)
ifdef HP3D_LIB_EXTRA_INCS
CCINCS   += $(HP3D_LIB_EXTRA_INCS)
endif

FC_WORK   = $(FC) $(CCDEFS) $(HP3D_MOD_OPTION_INC) $(CCINCS) \
            $(FFLAGS) $(EXTRA_FFLAGS)
CC_WORK   = $(CC) $(CDEFS) $(CCDEFS) $(CCINCS) $(CFLAGS) $(EXTRA_CFLAGS)

#-------------------------------------------------
PROB_FILE =

#-------------------------------------------------
# - HP3D make variable names
MK_COMMON_SRC   :=
MK_GMP_SRC      :=
MK_HP3D_SRC     :=

MK_COMMON_OBJS  :=
MK_GMP_OBJS     :=
MK_HP3D_OBJS    :=

-include ./sources.mk

# Convert source file paths to object file paths by replaying the base source
# directory with the object directory replacing suffix with .o
MK_COMMON_OBJS  += $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.c.o, \
                              $(filter %.c, $(MK_COMMON_SRC)))
MK_COMMON_OBJS  += $(patsubst $(SRC_PATH)/%.F, $(OBJ_PATH)/%.f.o, \
                              $(filter %.F, $(MK_COMMON_SRC)))
MK_COMMON_OBJS  += $(patsubst $(SRC_PATH)/%.F90, $(OBJ_PATH)/%.f90.o, \
                              $(filter %.F90, $(MK_COMMON_SRC)))

MK_GMP_OBJS     += $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.c.o, \
                              $(filter %.c, $(MK_GMP_SRC)))
MK_GMP_OBJS     += $(patsubst $(SRC_PATH)/%.F, $(OBJ_PATH)/%.f.o, \
                              $(filter %.F, $(MK_GMP_SRC)))
MK_GMP_OBJS     += $(patsubst $(SRC_PATH)/%.F90, $(OBJ_PATH)/%.f90.o, \
                              $(filter %.F90, $(MK_GMP_SRC)))

MK_HP3D_OBJS    += $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.c.o, \
                              $(filter %.c, $(MK_HP3D_SRC)))
MK_HP3D_OBJS    += $(patsubst $(SRC_PATH)/%.F, $(OBJ_PATH)/%.f.o, \
                              $(filter %.F, $(MK_HP3D_SRC)))
MK_HP3D_OBJS    += $(patsubst $(SRC_PATH)/%.F90, $(OBJ_PATH)/%.f90.o, \
                              $(filter %.F90, $(MK_HP3D_SRC)))

MK_ALL_HP3D_SRC  := $(MK_COMMON_SRC) \
		    $(MK_GMP_SRC) \
		    $(MK_HP3D_SRC)

MK_ALL_HP3D_OBJS := $(MK_COMMON_OBJS) \
		    $(MK_GMP_OBJS) \
		    $(MK_HP3D_OBJS)

# Introduce ranlib flag to suppress warnings on Mac OS
ARFLAG :=
ifeq ($(shell uname -s),Darwin)
	ARFLAG += -no_warning_for_no_symbols
endif

#-------------------------------------------------
# - Make -

.PHONY : lib

lib :
	@$(MAKE) depend;$(MAKE) depend;$(MAKE) fastlib;$(MAKE) install

#-------------------------------------------------
.PHONY : fast

fast :
	@$(MAKE) fastlib;$(MAKE) install

#-------------------------------------------------
.PHONY : fastlib

fastlib :
	@echo " - Building HP3D Library " $(LIB_PATH)/$(HP3D_LIBNAME)
	@echo "------------------------------------------"
	@$(MAKE) $(HP3D_LIBNAME)
	@echo "------------------------------------------"
	@echo " - Building Frontal Solver " $(FRONTAL_PATH)/$(FRONTAL_NAME)
	@(cd $(FRONTAL_PATH) ; $(MAKE) lib); \
	(cp -f $(FRONTAL_PATH)/$(FRONTAL_NAME) $(LIB_PATH));
	@echo "------------------------------------------"
	@echo " - COMPLEX, DEBUG         :" $(COMPLEX) $(DEBUG)
	@echo " - HP3D_USE_OPENMP        :" $(HP3D_USE_OPENMP)
	@echo " - HP3D_USE_INTEL_MKL     :" $(HP3D_USE_INTEL_MKL)
	@echo " - Preprocessor options   :" $(CCDEFS)
	@echo " - HP3D Library           :" $(LIB_PATH)/$(HP3D_LIBNAME)
	@echo " - Frontal Solver Library :" $(LIB_PATH)/$(FRONTAL_NAME)
	@echo " - Finished HP3D Library"
	@echo "------------------------------------------"
	@echo "    "

#-------------------------------------------------
$(HP3D_LIBNAME) : all
	@echo "Creating directory "$(LIB_PATH)
	@mkdir -p $(LIB_PATH)
	@echo "Creating library   "$(LIB_PATH)/$@
	@ar crS $(LIB_PATH)/$@ $(MK_ALL_HP3D_OBJS)
	@ranlib $(LIB_PATH)/$@ $(ARFLAG)

all : $(OBJ_PATH)/.dummy $(MK_ALL_HP3D_OBJS)

#-------------------------------------------------
$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_PATH); touch $@; \
	fi

-include ./depend_mod.mk

$(OBJ_PATH)/%.c.o : $(SRC_PATH)/%.c
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(CC_WORK) -o $@ -c $<

$(OBJ_PATH)/%.f.o : $(SRC_PATH)/%.F
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

$(OBJ_PATH)/%.f90.o : $(SRC_PATH)/%.F90
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

#-------------------------------------------------

source:
	@./build/sfmakedepend_src.sh $(SRC_PATH) sources.mk

depend: source
	@echo "Creating sources.mk and depend_mod.mk"
	@./build/sfmakedepend_mod.pl --silent $(subst -I,-I ,$(CCINCS)) \
	--moddir $(MODULE_PATH) \
	--objdir $(OBJ_PATH) \
	--srcprefix $(SRC_PATH) \
	--depend=mod \
	--longpath \
	--file depend_mod.mk \
	[file $(MK_ALL_HP3D_SRC) ]

install:
	@rm -rf $(INSTALL_PATH)
	@mkdir -p $(INSTALL_PATH)
	@mkdir -p $(INSTALL_PATH)/common
	@cp -r $(SRC_PATH)/common/* $(INSTALL_PATH)/common
	@mkdir -p $(INSTALL_PATH)/lib
	@if [ -f $(LIB_PATH)/$(HP3D_COMNAME) ]; then \
		cp $(LIB_PATH)/$(HP3D_COMNAME) $(INSTALL_PATH)/lib; \
	fi
	@if [ -f $(LIB_PATH)/$(HP3D_GMPNAME) ]; then \
		cp $(LIB_PATH)/$(HP3D_GMPNAME) $(INSTALL_PATH)/lib; \
	fi
	@if [ -f $(LIB_PATH)/$(HP3D_LIBNAME) ]; then \
		cp $(LIB_PATH)/$(HP3D_LIBNAME) $(INSTALL_PATH)/lib; \
	fi
	@if [ -f $(LIB_PATH)/$(FRONTAL_NAME) ]; then \
		cp $(LIB_PATH)/$(FRONTAL_NAME) $(INSTALL_PATH)/lib; \
	fi
	@mkdir -p $(INSTALL_PATH)/module;
	@cp $(MODULE_PATH)/*.mod $(INSTALL_PATH)/module

print:
	@echo $(MK_ALL_HP3D_SRC)

.PHONY : clean

clean:
	@echo "Clean frontal solver, object files, libraries, modules and mk files"
	@( cd ./$(SRC_PATH)/solver/frontal/base/real ; $(MAKE) clean )
	@( cd ./$(SRC_PATH)/solver/frontal/base/complex ; $(MAKE) clean )
	@rm -rf $(OBJ_PATH)
	@rm -rf $(OBJ_PATH_COMPLX) $(OBJ_PATH_REAL)
	@rm -rf $(HP3D_PATH)
	@rm -f $(LIB_PATH)/*.a
	@rm -f *~
	@rm -f $(MODULE_PATH)/*.mod
	@rm -f *.mk *.old *.mod
