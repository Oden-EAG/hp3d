FF=/opt/local/bin/gfortran-mp-4.7
F90=$(FF)
F77FLAGS= -O3
F90FLAGS= $(F77FLAGS)
LIB_PATH= -L/opt/local/lib
LIB_FFT= -ldrfftw -ldfftw
LIB_LAPACK= -llapack -lcblas -lf77blas -latlas


LIBS = $(LIB_PATH) $(LIB_FFT) $(LIB_LAPACK)

F77SR= xerror.f d1mach.f i1mach.f   zabs.f zacai.f zmlri.f zairy.f zrati.f zuni2.f zuni1.f zunk1.f zunk2.f zbunk.f zuchk.f zexp.f zdiv.f zsqrt.f zunhj.f zwrsk.f zseri.f zasyi.f zbuni.f zbinu.f zlog.f  zunik.f zuoik.f zmlt.f zs1s2.f zacon.f dgamln.f zkscl.f zshch.f  zbknu.f zbesh.f zbesj.f zbesy.f fdump.f j4save.f xgetua.f xercnt.f xerhlt.f xerprn.f xersve.f xermsg.f dxadd.f dxadj.f dxpsi.f dxpqnu.f dxpmu.f dxpmup.f dxpnrm.f dxqmu.f dxqnu.f dxred.f dxset.f dxlegf.f

F90SR=  fields-nmie.f90

MOM=$(F90SR:.f90=.o) $(F77SR:.f=.o)

.SUFFIXES:      .o .f .c .C .h .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<
.f.o:
	$(FF) $(F77FLAGS) -c $<

all:  test-code
test-code: $(MOM)
	$(F90) $(MOM) $(LIBS) -o $@

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.time
	rm -f *.mod
	rm -f core
	rm -f test-code
	rm -f *~
