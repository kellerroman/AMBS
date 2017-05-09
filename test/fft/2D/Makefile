FC=gfortran

LIBS = -lfftw3

INCLUDE = -I/usr/include/ 
#CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries #-Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

.PHONY: all compile run clean show

	
all: run show

compile: fftw_test_2d

run: compile 
	@./fftw_test_2d
show:
	@gnuplot -p plot.gnuplot
	
fftw_test_2d: fftw_test_2d.F90
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

clean:
	@rm -vf bin *.o *.mod *.h5 *.dat *.plt *.png *.xdmf *.txt *.bin fftw_test_2d *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~
