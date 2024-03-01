#------------------------------------------------
include ../../../../../m_options

SRC_PATH = .
OBJ_PATH = _obj_

CCINCS   = -I../../../../common/hp3d  \
	        -I../../../../../module

LIBS     =

FC_WORK  =  $(FC) $(CCINCS) $(FFLAGS)

SRC_FILE = bckwrdR.F90 compltR.F90 destR.F90 desvecR.F90 elmrhsR.F90 elmsolR.F90 \
	   frwcpR.F90 frwrsR.F90 loclsR.F90 locluR.F90 locrR.F90 resolR.F90 \
	   semrhsR.F90 surfspR.F90 surfssR.F90 symasmR.F90 symelmR.F90 \
	   unsasmR.F90 unselmR.F90 zdirioR.F90 zecondR.F90 zerorR.F90 \
	   zerosR.F90 zerouR.F90

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
	$(FC_WORK) -fPIC -o $@ -c $<

.PHONY : clean

clean :
	@rm -f $(OBJ_PATH)/*.o
	@rm -f *.o
	@rm -f *genmod.mod
	@rm -f *genmod.f90
	@rm -f libfrontal.a
