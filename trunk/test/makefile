info :
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

include ../m_options

COMPLEX         = 1
PARALLEL        = 0

CCDEFS          = -D C_MODE=1 \
                  -D PARALLEL_MODE=$(PARALLEL)

CCINCS          = -I../src/common/hp3d \
                  -I../src/common/graphics \
                  -I$(HP3D_MODULE) \
                  -I$(MUMPS_MPI_INC) \
                  -I./$(MODULE_PATH)


TEST  = random
TESTS = jacobi parameters physics random

HP3D_LIB     = ../complex/libhp3d.a ../complex/libhp3d_gmp.a ../complex/libhp3d_common.a

#LIBS         = $(HP3D_LIB)
#LIBS 	     =
FC_WORK      = $(FC)  $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

%.o : %.F
	$(FC_WORK) $(HP3D_COMMON) $(HP3D_MODULE) -o $@ -c $<

one : $(TEST).o
	$(FC_WORK) -o $(TEST) $(TEST).o  $(LDFLAGS) $(LIBS)

all :
	for e in $(TESTS) ; do \
		make one run TEST=$$e ;\
	done

run :
	./$(TEST)

clean :
	@/bin/rm -f *.o *~
	@for t in $(TESTS) ; do rm -f $$t ; done
