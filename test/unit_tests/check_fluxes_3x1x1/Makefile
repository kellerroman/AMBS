FC=gfortran

HOMEDIR=../../../

CGNSlibs = /usr/lib/libcgns.so -lhdf5
LIBS = $(CGNSlibs) -lm

INCLUDE = -I/usr/include/
CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries -Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

PRE_CONFIG=$(HOMEDIR)tools/bin/pre_config
PROGRAM=$(HOMEDIR)bin/AMBS

.PHONY: all cmp run clean cp grid_inf show

	
all: cmp run #show

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

	

run: data_in.cgns 
	@echo 'Running Grid-Adaper on simple Chamber'
	@$(PRE_CONFIG)
	@$(PROGRAM) 
show:
	@gnuplot -p plot.gnuplot

data_in.cgns: gridgen
	@./$<
	
gridgen: gridgen.F90
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)


clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~
