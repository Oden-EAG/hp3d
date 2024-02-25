#-------------------------------------------------
define HP3D_HELP_MSG

HP3D makefile targets:

  make
    Build and install the library
  make -j
    Build and install the library in parallel
  make clean
    Clean object/library files, problems and test directory
  make check
    Compile and run one test problem
  make test
    Compile and run test suite
  make problems
    Compile examples in problems directory
  make all
    Build/install the library, compile problems, and run tests
  make clean-lib
    Clean object and library files
  make clean-test
    Clean test directory
  make clean-problems
    Clean problems directory
  make clean-all
    Clean object/library files, problems and test directory

  Advanced usage:
    make lib
    make fast
    make fastlib
    make source
    make depend
    make install
    make print

endef

#-------------------------------------------------
# - Library path -

-include ./m_options

# Set variables if not set by m_options
HP3D_COMPLEX ?= 0
HP3D_DEBUG   ?= 0

HP3D_USE_INTEL_MKL ?= 0
HP3D_USE_MPI_F08   ?= 0
HP3D_USE_OPENMP    ?= 0
HP3D_USE_X11       ?= 1

####################################################

HP3D_LIBNAME   = libhp3d.a
FRONTAL_NAME   = libfrontal.a
HP3D_COMNAME   = libhp3d_common.a
HP3D_GMPNAME   = libhp3d_gmp.a

#-------------------------------------------------
# Complex version of library
ifeq ($(HP3D_COMPLEX),1)
INSTALL_PATH   = ./complex
CCINCS         =
FRONTAL_PATH   = ./$(SRC_PATH)/solver/frontal/base/complex
endif
#-------------------------------------------------
# Real version of library
ifeq ($(HP3D_COMPLEX),0)
INSTALL_PATH   = ./real
CCINCS         =
FRONTAL_PATH   = ./$(SRC_PATH)/solver/frontal/base/real
endif
#-------------------------------------------------
# - Compiler -
# -- Preprocessor options -- 
# Files with ending in .F and .F90 (not .f and .f90)
# are preprocessed according to certain directives
CCDEFS    = -D"HP3D_COMPLEX=$(HP3D_COMPLEX)" \
            -D"HP3D_DEBUG=$(HP3D_DEBUG)" \
            -D"HP3D_USE_INTEL_MKL=$(HP3D_USE_INTEL_MKL)" \
            -D"HP3D_USE_MPI_F08=$(HP3D_USE_MPI_F08)" \
            -D"HP3D_USE_OPENMP=$(HP3D_USE_OPENMP)" \
            -D"HP3D_USE_X11=$(HP3D_USE_X11)"

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
CCDEFS += -D"HP3D_USE_UNIX_DARWIN=1"
endif

#-------------------------------------------------
# - Make -
.PHONY : lib

lib :
	@$(MAKE) depend; $(MAKE) depend; $(MAKE) fastlib; $(MAKE) install

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
	@echo " - HP3D_COMPLEX           :" $(HP3D_COMPLEX)
	@echo " - HP3D_DEBUG             :" $(HP3D_DEBUG)
	@echo " - HP3D_USE_INTEL_MKL     :" $(HP3D_USE_INTEL_MKL)
	@echo " - HP3D_USE_MPI_F08       :" $(HP3D_USE_MPI_F08)
	@echo " - HP3D_USE_OPENMP        :" $(HP3D_USE_OPENMP)
	@echo " - HP3D_USE_X11           :" $(HP3D_USE_X11)
	@echo " - Preprocessor options   :" $(CCDEFS)
	@echo " - HP3D Library           :" $(LIB_PATH)/$(HP3D_LIBNAME)
	@echo " - Frontal Solver Library :" $(LIB_PATH)/$(FRONTAL_NAME)
	@echo " - Finished HP3D Library"
	@echo "------------------------------------------"
	@echo " "

#-------------------------------------------------
.PHONY : problems

problems: lib
	@echo "Compiling examples in 'problems' directory."
ifeq ($(HP3D_COMPLEX),1)
	@ERR=0; \
	echo "Compiling Maxwell Galerkin..."; \
	if ! $(MAKE) -C problems/MAXWELL/GALERKIN; then ERR=1; fi; \
	echo "Compiling Maxwell DPG Ultraweak..."; \
	if ! $(MAKE) -C problems/MAXWELL/ULTRAWEAK_DPG; then ERR=1; fi; \
	echo "Compiling Laser application..."; \
	if ! $(MAKE) -C problems/LASER/UW_COUPLED; then ERR=1; fi; \
	echo "Finished building complex-valued problems."; \
	if [ 0 -eq $${ERR} ]; then \
		echo "'make problems' passed."; \
	else \
		echo "'make problems' failed."; \
		exit 1; \
	fi ;
else
	@ERR=0; \
	echo "Compiling Poisson Galerkin..."; \
	if ! $(MAKE) -C problems/POISSON/GALERKIN; then ERR=1; fi; \
	echo "Compiling Poisson DPG Primal..."; \
	if ! $(MAKE) -C problems/POISSON/PRIMAL_DPG; then ERR=1; fi; \
	echo "Compiling Poisson DPG Ultraweak..."; \
	if ! $(MAKE) -C problems/POISSON/ULTRAWEAK_DPG; then ERR=1; fi; \
	echo "Compiling Vector Poisson Galerkin..."; \
	if ! $(MAKE) -C problems/VECTOR_POISSON/GALERKIN; then ERR=1; fi; \
	echo "Compiling Sheathed Hose application..."; \
	if ! $(MAKE) -C problems/SHEATHED_HOSE; then ERR=1; fi; \
	echo "Finished building real-valued problems."; \
	if [ 0 -eq $${ERR} ]; then \
		echo "'make problems' passed."; \
	else \
		echo "'make problems' failed."; \
		exit 1; \
	fi ;
endif

#-------------------------------------------------
.PHONY : check

check: lib
	@echo "Checking the hp3D library."
	@ERR=0; \
	if ! $(MAKE) -C test one run; then ERR=1; fi; \
	if ! $(MAKE) -C problems/MPI_TEST/TEST; then ERR=1; fi; \
	if [ 0 -eq $${ERR} ]; then \
		echo "'make check' passed."; \
	else \
		echo "'make check' failed."; \
		exit 1; \
	fi ;

#-------------------------------------------------
.PHONY : test

test: lib
	@echo "Testing the hp3D library."
	@if $(MAKE) -C test all; then \
		echo "'make test' passed."; \
	else \
		echo "'make test' failed."; \
		exit 1; \
	fi;

#-------------------------------------------------
.PHONY : all

all: test problems

#-------------------------------------------------
.PHONY : dummy

$(HP3D_LIBNAME) : dummy
	@echo "Creating directory "$(LIB_PATH)
	@mkdir -p $(LIB_PATH)
	@echo "Creating library   "$(LIB_PATH)/$@
	@ar crS $(LIB_PATH)/$@ $(MK_ALL_HP3D_OBJS)
	@ranlib $(LIB_PATH)/$@ $(ARFLAG)

dummy : $(OBJ_PATH)/.dummy $(MK_ALL_HP3D_OBJS)

#-------------------------------------------------
$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@mkdir  -p $(OBJ_PATH); touch $@;

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
.PHONY : source depend install

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

#-------------------------------------------------
.PHONY : help

help:
	$(info $(value HP3D_HELP_MSG))
	@true

#-------------------------------------------------
.PHONY : print

print:
	@echo $(MK_ALL_HP3D_SRC)

#-------------------------------------------------
.PHONY : clean

clean: clean-all

#-------------------------------------------------
.PHONY : clean-lib

clean-lib:
	@echo "Clean frontal solver, object files, libraries, modules and mk files"
	@( $(MAKE) -C ./$(SRC_PATH)/solver/frontal/base/real    clean )
	@( $(MAKE) -C ./$(SRC_PATH)/solver/frontal/base/complex clean )
	@rm -rf $(OBJ_PATH_COMPLEX) $(OBJ_PATH_REAL)
	@rm -rf $(HP3D_PATH_COMPLEX) $(HP3D_PATH_REAL)
	@rm -f $(LIB_PATH)/*.a
	@rm -f *~
	@rm -f $(MODULE_PATH)/*.mod
	@rm -f $(MODULE_PATH)/*genmod.f90
	@rm -f *.mk *.old *.mod

#-------------------------------------------------
.PHONY : clean-test

clean-test:
	@echo "Clean 'test' directory."
	@($(MAKE) -C test clean);

#-------------------------------------------------
.PHONY : clean-problems

clean-problems:
	@echo "Clean 'problems' directory."
	@($(MAKE) -C problems/MPI_TEST/TEST clean);
ifeq ($(HP3D_COMPLEX),1)
	@($(MAKE) -C problems/MAXWELL/GALERKIN      clean);
	@($(MAKE) -C problems/MAXWELL/ULTRAWEAK_DPG clean);
	@($(MAKE) -C problems/LASER/UW_COUPLED      clean);
else
	@($(MAKE) -C problems/POISSON/GALERKIN        clean);
	@($(MAKE) -C problems/POISSON/PRIMAL_DPG      clean);
	@($(MAKE) -C problems/POISSON/ULTRAWEAK_DPG   clean);
	@($(MAKE) -C problems/VECTOR_POISSON/GALERKIN clean);
	@($(MAKE) -C problems/SHEATHED_HOSE           clean);
endif

#-------------------------------------------------
.PHONY : clean-all

clean-all: clean-lib clean-test clean-problems
