info :
	@echo "make all"
	@echo "make one TEST=test"
	@echo "make run TEST=test"
	@echo "    test="
	@for t in $(TESTS); do echo "         "$$t; done

include ../m_options

OBJ_PATH     = _obj_
EXEC_PATH    = _exec_
MODULE_PATH  = $(OBJ_PATH)

SRC_FILE     = locate.F90           \
               hexa_adapt_h.F90     \
               parameters.F90       \
               random_refine.F90

OBJ          = $(addprefix $(OBJ_PATH)/, $(SRC_FILE:.F90=.o))

# Tests
TEST  = hexa_adapt_h
TESTS = locate          \
        hexa_adapt_h    \
        parameters      \
        random_refine

# Includes and defs
CCDEFS = $(PROB_PP_DEFS)
CCINCS = $(PROB_INCS)

# User defined libraries
LIBS   = $(PROB_LIBS)

# Compiler flags
FC_WORK = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

#-------------------------------------------------

$(OBJ) : $(OBJ_PATH)/%.o : %.F90
	@echo ""
	@echo "Compiling $<"
	$(FC_WORK) -o $@ -c $<
	
$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@mkdir -p $(OBJ_PATH); touch $@;
	@mkdir -p $(EXEC_PATH); touch $@;

#-------------------------------------------------
.PHONY : one
one : $(OBJ_PATH)/.dummy $(OBJ_PATH)/$(TEST).o
	@echo ""
	@echo "- Linking - " $(TEST)
	$(FC_WORK) -o $(EXEC_PATH)/$(TEST) $(OBJ_PATH)/$(TEST).o $(LIBS) $(LDFLAGS)

#-------------------------------------------------
.PHONY : all
all :
	@for e in $(TESTS) ; do \
		make one run TEST=$$e ;\
	done

#-------------------------------------------------
.PHONY : run
run : one
	@echo ""
	@echo "Running test: "$(TEST)
	@./$(EXEC_PATH)/$(TEST)

#-------------------------------------------------
.PHONY : clean
clean :
	@/bin/rm -f *.o *~
	@/bin/rm -rf $(OBJ_PATH)
	@/bin/rm -rf $(EXEC_PATH)
	@for t in $(TESTS) ; do rm -f $$t ; done
