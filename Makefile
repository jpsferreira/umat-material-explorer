# UMAT Material Explorer
# Explore material behavior under different loading scenarios
# and fit material parameters to experimental data.

.PHONY: all cyclic monotonic fitting clean

all: cyclic monotonic fitting

cyclic:
	$(MAKE) -C cyclic

monotonic:
	$(MAKE) -C monotonic

fitting:
	$(MAKE) -C fitting

clean:
	$(MAKE) -C cyclic clean
	$(MAKE) -C monotonic clean
	$(MAKE) -C fitting clean
