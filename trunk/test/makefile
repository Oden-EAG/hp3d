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

SRC_FILE     = cross_product.F90         \
               decode.F90                \
               dot_product.F90           \
               encod_decod.F90           \
               encode_decode_order.F90   \
               encodg_decodg.F90         \
               fact.F90                  \
               hexa_adapt_h_aniso.F90    \
               hexa_adapt_h_case1.F90    \
               hexa_adapt_h_iso.F90      \
               invert.F90                \
               isgn.F90                  \
               locate_char.F90           \
               locate.F90                \
               near.F90                  \
               parameters.F90            \
               real_imag.F90

OBJ          = $(addprefix $(OBJ_PATH)/, $(SRC_FILE:.F90=.o))

# Tests
TEST     = decode
TEST_SEQ = cross_product             \
           decode                    \
           dot_product               \
           encod_decod               \
           encode_decode_order       \
           encodg_decodg             \
           fact                      \
           invert                    \
           isgn                      \
           locate_char               \
           locate                    \
           near                      \
           parameters                \
           real_imag
TEST_PAR = hexa_adapt_h_aniso        \
           hexa_adapt_h_case1        \
           hexa_adapt_h_iso

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
.PHONY : seq
seq :
	@for e in $(TEST_SEQ) ; do \
		make one run TEST=$$e ;\
	done

#-------------------------------------------------
.PHONY : par
par :
	@for e in $(TEST_PAR) ; do \
		make one run-par TEST=$$e ;\
	done

#-------------------------------------------------
.PHONY : all
all : seq par
	@echo ""

#-------------------------------------------------
.PHONY : run
run : one
	@echo ""
	@echo "Running test: "$(TEST)
	@./$(EXEC_PATH)/$(TEST)

#-------------------------------------------------
.PHONY : run-par
run-par : one
	@echo ""
	@echo "Running test: "$(TEST)
	@mpirun -np 4 ./$(EXEC_PATH)/$(TEST)

#-------------------------------------------------
.PHONY : clean
clean :
	@/bin/rm -f *.o *~
	@/bin/rm -rf $(OBJ_PATH)
	@/bin/rm -rf $(EXEC_PATH)
	@for t in $(TEST_SEQ) ; do rm -f $$t ; done
	@for t in $(TEST_PAR) ; do rm -f $$t ; done
