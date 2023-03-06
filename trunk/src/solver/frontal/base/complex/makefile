#------------------------------------------------
include ../../../../../m_options

SRC_PATH = .
OBJ_PATH = _obj_

CCINCS   = -I../../../../common/hp3d  \
	   -I../../../../common/frontal/complex

LIBS     =

FC_WORK  =  $(FC) $(CCINCS) $(FFLAGS)

SRC_FILE = bckwrd.F complt.F dest.F desvec.F elmrhs.F elmsol.F \
	   frwcp.F frwrs.F frwsr.F locls.F loclu.F locr.F \
	   resol.F semrhs.F surfsp.F surfss.F symasm.F \
	   symelm.F unsasm.F unselm.F zdirio.F zdiriodum.F \
	   zecond.F zeror.F zeros.F zerou.F

SRC      = $(addprefix $(SRC_PATH)/.$(SRC_FILE))
OBJ      = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F=.o))

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
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_PATH); touch $@; \
	fi

$(OBJ_PATH)/%.o : $(SRC_PATH)/%.F
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<

.PHONY : clean

clean :
	@rm -f $(OBJ_PATH)/*.o
	@rm -f *.o
	@rm -f libfrontal.a
