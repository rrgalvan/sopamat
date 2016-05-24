PGFORTRAN_DIR = /home/rrgalvan/local/linux86-64/2016/bin
FC = $(PGFORTRAN_DIR)/pgfortran
FCFLAGS = -Minform=warn -acc -Minfo -ta=nvidia,time

TARGETS = finite_diff_1d

all: $(TARGETS)

finite_diff_1d: finite_diff_1d.f90
	$(FC) $(FCFLAGS) $< -o $@

clean:
	rm -f $(TARGETS) *~

run:
	./finite_diff_1d
