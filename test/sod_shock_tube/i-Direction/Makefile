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

str.dx: ../stubex
	@./$<
	
../stubex: ../stubex.f
	@$(FC) $(CFLAGS) $< -o $@

run: data_in.cgns str.dx
	@echo 'Running Grid-Adaper on simple Chamber'
	@$(PRE_CONFIG)
	@$(PROGRAM) #-debug
	@$(HOMEDIR)/tools/bin/extract_1D 
show:
	@gnuplot -p plot.gnuplot

data_in.cgns: gridgen
	@./$< 101
	
gridgen: gridgen.F90 cgnslib.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS)
	
cgnslib.o: $(HOMEDIR)src/cgns.F90
	$(FC) -c $(CFLAGS) $< -o $@ 

grid_rot: cmp gridgen
	./gridgen 101 0.0
	$(MAKE) run
	cp data_1d.csv data_1d_00.csv
	./gridgen 101 45.0
	$(MAKE) run
	cp data_1d.csv data_1d_45.csv
	./gridgen 101 90.0
	$(MAKE) run
	cp data_1d.csv data_1d_90.csv
	@gnuplot -p plot_grid_rot.gnuplot

clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
