FC=gfortran

HOMEDIR=../../../

CGNSlibs = /usr/lib/libcgns.so -lhdf5
LIBS = $(CGNSlibs) -lm
CFLAGS += -O3
CFLAGS += -fbounds-check -fcheck-array-temporaries -Wall
#CFLAGS += -DDEBUG
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
	
gridgen: gridgen.F90 cgnslib.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS)
	
cgnslib.o: $(HOMEDIR)src/cgns.F90
	$(FC) -c $(CFLAGS) $< -o $@ 

clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
