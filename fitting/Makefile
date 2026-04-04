FC      = gfortran
FFLAGS  = -O2 -Wall

SRCS    = umat.for uexternaldb.for uexternald.for sga.f95 main.f90
INCS    = param_umat.inc aba_param.inc params.inc

all: sga_fitting

sga_fitting: $(SRCS) $(INCS)
	$(FC) $(FFLAGS) -o $@ $(SRCS)

run: sga_fitting
	./sga_fitting

clean:
	rm -f sga_fitting *.o *.mod
