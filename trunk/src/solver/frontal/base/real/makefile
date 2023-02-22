#------------------------------------------------
include ../../../../../m_options

SRC_PATH = .
OBJ_PATH = _obj_

CCINCS   = -I../../../../common/hp3d  \
	   -I../../../../common/frontal/real

LIBS     =

FC_WORK  =  $(FC) $(CCINCS) $(FFLAGS)

SRC_FILE = bckwrdR.F compltR.F destR.F desvecR.F elmrhsR.F elmsolR.F \
	   frwcpR.F frwrsR.F loclsR.F locluR.F locrR.F resolR.F \
	   semrhsR.F surfspR.F surfssR.F symasmR.F symelmR.F \
	   unsasmR.F unselmR.F zdirioR.F zecondR.F zerorR.F \
	   zerosR.F zerouR.F

SRC      = $(addprefix $(SRC_PATH)/.$(SRC_FILE))
OBJ      = $(addprefix $(OBJ_PATH)/,$(SRC_FILE:.F=.o))

.PHONY : lib

lib : libfrontal.a
	@echo "Creating library "$@/$^

libfrontal.a : $(OBJ_PATH)/.dummy $(OBJ)
	@ar crS $@ $(OBJ)
	@ranlib $@ -no_warning_for_no_symbols

$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@if [ -d $(OBJ_PATH) ]; then \
		touch $@; \
	else mkdir $(OBJ_PATH); touch $@; \
	fi

$(OBJ_PATH)/%.o : $(SRC_PATH)/%.F
	@echo "Compiling $<"
	$(FC_WORK) -fPIC -o $@ -c $<

.PHONY : clean

clean :
	@rm -f $(OBJ_PATH)/*.o
	@rm -f *.o
	@rm -f libfrontal.a
