FC      = gfortran
FFLAGS  = -O2
UMAT    = ../umat/umat.for ../umat/resetdfgr.for

all: cyclic

cyclic: cyclic.f90 $(UMAT)
	$(FC) $(FFLAGS) -I../umat -o $@ $(UMAT) $<

run: cyclic
	ln -sf ../umat/fibers.inp fibers.inp
	./cyclic

clean:
	rm -f cyclic *.o *.mod fibers.inp
	rm -rf stress_curves
