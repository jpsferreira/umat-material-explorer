FC      = gfortran
FFLAGS  = -O2
UMAT    = ../umat/umat.for ../umat/resetdfgr.for

all: monotonic

monotonic: monotonic.f90 $(UMAT)
	$(FC) $(FFLAGS) -I../umat -o $@ $(UMAT) $<

run: monotonic
	ln -sf ../umat/fibers.inp fibers.inp
	./monotonic

clean:
	rm -f monotonic *.o *.mod fibers.inp
	rm -rf stress_curves
