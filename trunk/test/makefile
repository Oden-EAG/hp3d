info :
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

include ../m_options

COMPLEX         = 0
CCDEFS          = -D C_MODE=0
CCINCS          = -I../src/common/hp3d \
                  -I../src/common/graphics \
                  -I../$(MODULE_PATH)

TEST  = locate
TESTS = locate parameters random

HP3D_LIB     = ../real/lib/libhp3d.a
#HP3D_LIB     = ../complex/lib/libhp3d.a

LIBS         = $(HP3D_LIB)
#LIBS 	     =
FC_WORK      = $(FC) $(CCDEFS) $(CCINCS) $(FFLAGS) $(EXTRA_FFLAGS)

%.o : %.F
	$(FC_WORK) $(HP3D_LIB) -o $@ -c $<

one : $(TEST).o
	$(FC_WORK) -o $(TEST) $(TEST).o $(LIBS)

all :
	for e in $(TESTS) ; do \
		make one run TEST=$$e ;\
	done

run :
	./$(TEST)

clean :
	@/bin/rm -f *.o *~
	@for t in $(TESTS) ; do rm -f $$t ; done
