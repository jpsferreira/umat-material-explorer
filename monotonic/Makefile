FC      = gfortran
FFLAGS  = -O2 -Wall
UMAT    = ../umat/umat.for ../umat/resetdfgr.for
INCS    = ../umat/param_umat.inc ../umat/aba_param.inc

all: monotonic

monotonic: monotonic.f90 $(UMAT) $(INCS)
	$(FC) $(FFLAGS) -I../umat -o $@ $(UMAT) $<

run: monotonic
	ln -sf ../umat/fibers.inp fibers.inp
	./monotonic

clean:
	rm -f monotonic *.o *.mod fibers.inp
	rm -rf stress_curves
