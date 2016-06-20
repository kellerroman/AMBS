FC=gfortran

HOMEDIR=../../

CGNSlibs = /usr/lib/libcgns.so -lhdf5
LIBS = $(CGNSlibs) -lm

INCLUDE = -I/usr/include/
CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries -Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

PRE_CONFIG=$(HOMEDIR)tools/bin/pre_config
EXTRACT_1D=$(HOMEDIR)tools/bin/extract_1D
PROGRAM=$(HOMEDIR)bin/AMBS

.PHONY: all cmp run clean cp grid_inf show

	
all: cmp run #show

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

	

run: data_in.cgns config.bin
	@echo 'Running Grid-Adaper on simple Chamber'
	@$(PROGRAM) #-debug
	@$(EXTRACT_1D)
show:
	@gnuplot -p plot.gnuplot

data_in.cgns: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)
gridgen: gridgen.F90 cgnslib.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)
	
cgnslib.o: $(HOMEDIR)src/cgns.F90
	$(FC) -c $(CFLAGS) $< -o $@ 


clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
