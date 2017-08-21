HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

CALC_ERROR=../calc_error &>> errors.txt

.PHONY: all run clean grid_inf show
.DEFAULT_GOAL := all

	
all: compile_global run #show

str.dx: ../stubex
	@./$<
	
../stubex: ../stubex.f
	@$(FC) $(CFLAGS) $< -o $@

run: data_in.h5 str.dx config.bin
	@echo 'Running Grid-Adaper on simple Chamber'
	@$(PROGRAM) #-debug
	@$(EXTRACT_1D)
	@$(CREATE_XDMF)

show:
	@gnuplot -p plot.gnuplot

data_in.h5: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)

gridgen: gridgen.F90 $(HOMEDIR)tools/obj/module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

sens: compile_global gridgen str.dx
	./gridgen 101
	@$(PRE_CONFIG) disc/space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_100_sp_1.csv
	@$(PRE_CONFIG) disc/space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_100_sp_2.csv
	./gridgen 201
	@$(PRE_CONFIG) disc/space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_200_sp_1.csv
	@$(PRE_CONFIG) disc/space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_200_sp_2.csv
	./gridgen 501
	@$(PRE_CONFIG) disc/space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_500_sp_1.csv
	@$(PRE_CONFIG) disc/space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@mv data_1d.csv data_ni_500_sp_2.csv
	@gnuplot -p plot_sens.gnuplot

grid_rot: compile_global gridgen
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
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex *.h5 *.csv *.xdmf
	@rm -vrf *~
