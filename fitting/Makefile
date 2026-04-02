FC      = gfortran
FFLAGS  = -O2

SRCS    = umat_gho.for uexternaldb.for uexternald.for sga.f95 main.f90

all: sga_fitting

sga_fitting: $(SRCS)
	$(FC) $(FFLAGS) -o $@ $^

run: sga_fitting
	./sga_fitting

clean:
	rm -f sga_fitting *.o *.mod
