FC      = gfortran
FFLAGS  = -O2 -Wall
UMAT    = ../umat/umat.for ../umat/resetdfgr.for
INCS    = ../umat/param_umat.inc ../umat/aba_param.inc

all: cyclic

cyclic: cyclic.f90 $(UMAT) $(INCS)
	$(FC) $(FFLAGS) -I../umat -o $@ $(UMAT) $<

run: cyclic
	ln -sf ../umat/fibers.inp fibers.inp
	./cyclic

clean:
	rm -f cyclic *.o *.mod fibers.inp
	rm -rf stress_curves
