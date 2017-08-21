HOMEDIR=../../

include $(HOMEDIR)Makefile.global

INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

.PHONY: all run clean grid_rot show
.DEFAULT_GOAL := all

	
all: compile_global run #show

run: data_in.h5 config.bin
	@echo ' Running Shu Osher Shock Problem'
	@$(PROGRAM) #-debug
	@$(CREATE_XDMF)
	@$(EXTRACT_1D)
show:
	@gnuplot -p plot.gnuplot

data_in.h5: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

sens: cmp gridgen 
	@rm -f errors.txt
	./gridgen 1001
	@$(PRE_CONFIG) disc/space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_100_sp_1.csv
	@$(PRE_CONFIG) disc/space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_100_sp_2.csv
	./gridgen 2001
	@$(PRE_CONFIG) disc/space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_200_sp_1.csv
	@$(PRE_CONFIG) disc/space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_200_sp_2.csv
	./gridgen 5001
	@$(PRE_CONFIG) disc/space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_500_sp_1.csv
	@$(PRE_CONFIG) disc/space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_500_sp_2.csv
	./gridgen 10001
	@$(PRE_CONFIG) disc/space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_1000_sp_1.csv
	@$(PRE_CONFIG) disc/space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_1000_sp_2.csv
	@gnuplot -p plot_sens.gnuplot
        
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
	@rm -vf *.o *.mod *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.dx stubex *.csv *.h5 *.xdmf
	@rm -vrf *~
