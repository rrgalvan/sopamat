PGFORTRAN_DIR = /home/rrgalvan/local/linux86-64/2016/bin
FC = $(PGFORTRAN_DIR)/pgfortran
FCFLAGS = -Minform=warn -acc -Minfo

TARGETS = openacc_ex1 openacc_ex2

all: $(TARGETS)

openacc_ex1: openacc_ex1.f90
	$(FC) $(FCFLAGS) $< -o $@

openacc_ex2: openacc_ex2.f90
	$(FC) $(FCFLAGS) $< -o $@

clean:
	rm -f $(TARGETS) *~
