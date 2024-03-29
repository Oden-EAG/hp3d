info :
	@echo "make all"
	@echo "make one TEST=test"
	@echo "make run TEST=test"
	@echo "    test="
	@for t in $(TESTS); do echo "         "$$t; done

include ../m_options

HP3D_MPIEXEC    ?= mpirun
HP3D_MPIEXEC_NP ?= -np
HP3D_MPI_NP     ?= 2

HP3D_MPI_SEQ ?= $(HP3D_MPIEXEC) $(HP3D_MPIEXEC_NP) 1
HP3D_MPI_PAR ?= $(HP3D_MPIEXEC) $(HP3D_MPIEXEC_NP) $(HP3D_MPI_NP)

OBJ_PATH     = _obj_
EXEC_PATH    = _exec_
MODULE_PATH  = $(OBJ_PATH)

SRC_FILE     = attr_to_index.F90         \
               cross_product.F90         \
               decode.F90                \
               dot_product.F90           \
               encod_decod.F90           \
               encode_decode_order.F90   \
               encodg_decodg.F90         \
               fact.F90                  \
               hexa_adapt_h_aniso.F90    \
               hexa_adapt_h_case1.F90    \
               hexa_adapt_h_iso.F90      \
               ij_to_packed.F90          \
               index_to_attr.F90         \
               invert.F90                \
               isgn.F90                  \
               locate_char.F90           \
               locate.F90                \
               near.F90                  \
               parameters.F90            \
               real_imag.F90
ifeq ($(HP3D_COMPLEX),1)
SRC_FILE    += conv_maxw.F90             \
               dump_maxw.F90             \
               poly_maxw.F90
else
SRC_FILE    += conv_pois.F90             \
               dump_pois.F90             \
               poly_pois.F90
endif

OBJ          = $(addprefix $(OBJ_PATH)/, $(SRC_FILE:.F90=.o))

# Tests
TEST     = decode
TEST_SEQ = attr_to_index             \
           cross_product             \
           decode                    \
           dot_product               \
           encod_decod               \
           encode_decode_order       \
           encodg_decodg             \
           fact                      \
           ij_to_packed              \
           index_to_attr             \
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
ifeq ($(HP3D_COMPLEX),1)
TEST_SEQ+= conv_maxw                 \
           dump_maxw                 \
           poly_maxw
else
TEST_SEQ+= conv_pois                 \
           dump_pois                 \
           poly_pois
endif

# Includes and defs
CCDEFS = $(PROB_PP_DEFS)
CCINCS = $(PROB_INCS)

# User defined libraries
LIBS   = $(PROB_LIBS)

# Compiler flags
FC_WORK = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

#-------------------------------------------------

$(OBJ) : $(OBJ_PATH)/%.o : %.F90
	@printf "\nCompiling $<\n"
	$(FC_WORK) -o $@ -c $<
	
$(OBJ_PATH)/.dummy :
	@echo "Creating directory " $(OBJ_PATH)
	@mkdir -p $(OBJ_PATH); touch $@;
	@mkdir -p $(EXEC_PATH); touch $@;

#-------------------------------------------------
.PHONY : one
one : $(OBJ_PATH)/.dummy $(OBJ_PATH)/$(TEST).o
	@printf "\n - Linking - "$(TEST)"\n"
	$(FC_WORK) -o $(EXEC_PATH)/$(TEST) $(OBJ_PATH)/$(TEST).o $(LIBS) $(LDFLAGS)

#-------------------------------------------------
.PHONY : seq
seq :
	@seq_pass=0; seq_fail=0; \
	for e in $(TEST_SEQ) ; do \
		pass=1; $(MAKE) one run TEST=$$e || pass=0; \
		if [ $$((pass)) -eq 1 ]; then \
			seq_pass=$$((seq_pass+1)); \
		else \
			seq_fail=$$((seq_fail+1)); \
		fi ; \
	done ; \
	if [ $$((seq_fail)) -eq 0 ]; then \
		printf "\n All "$$((seq_pass))"%s\n\n" \
			    " Sequential Tests Passed."; \
	else \
		printf "\n     "$$((seq_pass))/$$((seq_pass+seq_fail))"%s\n" \
			    " Sequential Tests Passed."; \
		printf   "     "$$((seq_fail))/$$((seq_pass+seq_fail))"%s\n\n" \
			    " Sequential Tests FAILED!"; \
		exit 1; \
	fi ; \
	echo $$((seq_pass))" "$$((seq_fail)) > $(OBJ_PATH)/.seq_test

#-------------------------------------------------
.PHONY : par
par :
	@par_pass=0; par_fail=0; \
	for e in $(TEST_PAR) ; do \
		pass=1; $(MAKE) one run-par TEST=$$e || pass=0; \
		if [ $$((pass)) -eq 1 ]; then \
			par_pass=$$((par_pass+1)); \
		else \
			par_fail=$$((par_fail+1)); \
		fi ; \
	done ; \
	if [ $$((par_fail)) -eq 0 ]; then \
		printf "\n All "$$((par_pass))"%s\n\n" \
			    " Parallel Tests Passed."; \
	else \
		printf "\n     "$$((par_pass))/$$((par_pass+par_fail))"%s\n" \
			    " Parallel Tests Passed."; \
		printf   "     "$$((par_fail))/$$((par_pass+par_fail))"%s\n\n" \
			    " Parallel Tests FAILED!"; \
		exit 1; \
	fi ; \
	echo $$((par_pass))" "$$((par_fail)) > $(OBJ_PATH)/.par_test

#-------------------------------------------------
.PHONY : all
all : seq par
	@read -r seq_pass seq_fail < $(OBJ_PATH)/.seq_test; \
	 read -r par_pass par_fail < $(OBJ_PATH)/.par_test; \
	 pass=$$((seq_pass+par_pass)); \
	 fail=$$((seq_fail+par_fail)); \
	 echo "++++++++++++++++++++++++++++++++++++++++"; \
	if [ $$((fail)) -eq 0 ]; then \
		printf "\n All "$$((pass))"%s\n\n" \
			    " Tests Passed."; \
	else \
		printf "\n     "$$((pass))/$$((pass+fail))"%s\n" \
			    " Tests Passed."; \
		printf   "     "$$((fail))/$$((pass+fail))"%s\n\n" \
			    " Tests FAILED!"; \
		exit 1; \
	fi ; \
	 echo "++++++++++++++++++++++++++++++++++++++++";

#-------------------------------------------------
.PHONY : run
run : one
	@printf "\n Running test: "$(TEST)"\n"
	@pass=1; $(HP3D_MPI_SEQ) ./$(EXEC_PATH)/$(TEST) || pass=0; \
	if [ $$((pass)) -eq 1 ]; then \
		printf "          test_"$(TEST)" PASSED\n\n"; \
		true; \
	else \
		printf "          test_"$(TEST)" FAILED\n\n"; \
		false; \
		exit 1; \
	fi ;

#-------------------------------------------------
.PHONY : run-par
run-par : one
	@printf "\n Running test: "$(TEST)"\n"
	@pass=1; $(HP3D_MPI_PAR) ./$(EXEC_PATH)/$(TEST) || pass=0; \
	if [ $$((pass)) -eq 1 ]; then \
		printf "          test_"$(TEST)" PASSED\n\n"; \
		true; \
	else \
		printf "          test_"$(TEST)" FAILED\n\n"; \
		false; \
		exit 1; \
	fi ;

#-------------------------------------------------
.PHONY : clean
clean :
	@/bin/rm -f *.o *~
	@/bin/rm -rf $(OBJ_PATH)
	@/bin/rm -rf $(EXEC_PATH)
	@for t in $(TEST_SEQ) ; do rm -f $$t ; done
	@for t in $(TEST_PAR) ; do rm -f $$t ; done
