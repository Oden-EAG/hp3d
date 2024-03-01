#------------------------------------------------
include ../../../../../m_options

SRC_PATH = .
OBJ_PATH = _obj_

CCINCS   = -I../../../../common/hp3d  \
	        -I../../../../../module

LIBS     =

FC_WORK  =  $(FC) $(CCINCS) $(FFLAGS)

SRC_FILE = bckwrd.F90 complt.F90 dest.F90 desvec.F90 elmrhs.F90 elmsol.F90 \
	   frwcp.F90 frwrs.F90 locls.F90 loclu.F90 locr.F90 \
	   resol.F90 semrhs.F90 surfsp.F90 surfss.F90 symasm.F90 \
	   symelm.F90 unsasm.F90 unselm.F90 zdirio.F90 zdiriodum.F90 \
	   zecond.F90 zeror.F90 zeros.F90 zerou.F90

SRC      = $(addprefix $(SRC_PATH)/.$(SRC_FILE))
OBJ      = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F90=.o))

# Introduce ranlib flag to suppress warnings on Mac OS
ARFLAG :=
ifeq ($(shell uname -s),Darwin)
	ARFLAG += -no_warning_for_no_symbols
endif

.PHONY : lib

lib : libfrontal.a
	@echo "Creating library "$@/$^

libfrontal.a : $(OBJ_PATH)/.dummy $(OBJ)
	@ar crS $@ $(OBJ)
	@ranlib $@ $(ARFLAG)

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@mkdir -p $(OBJ_PATH); touch $@;

$(OBJ_PATH)/%.o : $(SRC_PATH)/%.F90
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

.PHONY : clean

clean :
	@rm -f $(OBJ_PATH)/*.o
	@rm -f *.o
	@rm -f *genmod.mod
	@rm -f *genmod.f90
	@rm -f libfrontal.a
