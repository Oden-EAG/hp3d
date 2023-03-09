info :
	@echo "make all"
	@echo "make one TEST=test"
	@echo "make run TEST=test"
	@echo "    test="
	@for t in $(TESTS); do echo "         "$$t; done

include ../m_options

OBJ_PATH     = _obj_
MODULE_PATH  = $(OBJ_PATH)

SRC_FILE     = locate.F \
               parameters.F \
               random_refine.F

OBJ          = $(addprefix $(OBJ_PATH)/, $(SRC_FILE:.F=.o))

# Tests
TEST  = locate
TESTS = locate parameters random_refine

# Includes and defs
CCDEFS = $(PROB_PP_DEFS)
CCINCS = $(PROB_INCS)

# User defined libraries
LIBS   = $(PROB_LIBS)

# Compiler flags
FC_WORK = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

#-------------------------------------------------

$(OBJ) : $(OBJ_PATH)/%.o : %.F
	@echo ""
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<
	
$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@mkdir -p $(OBJ_PATH); touch $@;

one : $(OBJ_PATH)/.dummy $(OBJ_PATH)/$(TEST).o
	@echo ""
	@echo "- Linking - " $(TEST)
	$(FC_WORK) -o $(TEST) $(OBJ_PATH)/$(TEST).o $(LIBS) $(LDFLAGS)
	@echo $(CCINCS)

all :
	for e in $(TESTS) ; do \
		make one run TEST=$$e ;\
	done

run : one
	@echo ""
	@echo "Running test: "$(TEST)
	./$(TEST)

clean :
	@/bin/rm -f *.o *~
	@/bin/rm -rf $(OBJ_PATH)
	@for t in $(TESTS) ; do rm -f $$t ; done
