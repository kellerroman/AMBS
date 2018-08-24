HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

CALC_ERROR=../calc_error &>> errors.txt

.PHONY: all run clean grid_inf show
.DEFAULT_GOAL := all

	
all: compile run #show

compile:compile_global
	@$(MAKE) -C ../ calc_error

str.dx: ../stubex
	@./$<
	
../stubex: ../stubex.f
	@$(FC) $(CFLAGS) $< -o $@

run: data_in.h5 str.dx config.bin
	@echo ' Running SOD Shock Tube'
	@$(PROGRAM) #-debug
	@$(CREATE_XDMF)
	@$(EXTRACT_1D)

show:
	@gnuplot -p plot.gnuplot

data_in.h5: gridgen
	@./$< 101

config.bin: config.cfg	
	@$(PRE_CONFIG)

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

lim: compile gridgen str.dx
	@rm -f errors.txt
	@./gridgen 501
	@$(PRE_CONFIG) space_order=2 limiter=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_lim=1.csv
	@mv error.csv error_lim=1.csv
	@$(PRE_CONFIG) space_order=2 limiter=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_lim=2.csv
	@mv error.csv error_lim=2.csv
	@$(PRE_CONFIG) space_order=2 limiter=3
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_lim=3.csv
	@mv error.csv error_lim=3.csv
	@gnuplot -p plot_lim.gnuplot
	@cat errors.txt

sens: compile gridgen str.dx
	@rm -f errors.txt
	@./gridgen 51
	@$(PRE_CONFIG) space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_050_sp_1.csv
	@mv error.csv error_ni_050_sp_1.csv
	@$(PRE_CONFIG) space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_050_sp_2.csv
	@mv error.csv error_ni_050_sp_2.csv
	@./gridgen 101
	@$(PRE_CONFIG) space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_100_sp_1.csv
	@mv error.csv error_ni_100_sp_1.csv
	@$(PRE_CONFIG) space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_100_sp_2.csv
	@mv error.csv error_ni_100_sp_2.csv
	@./gridgen 201
	@$(PRE_CONFIG) space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_200_sp_1.csv
	@mv error.csv error_ni_200_sp_1.csv
	@$(PRE_CONFIG) space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_200_sp_2.csv
	@mv error.csv error_ni_200_sp_2.csv
	@./gridgen 401
	@$(PRE_CONFIG) space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_400_sp_1.csv
	@mv error.csv error_ni_400_sp_1.csv
	@$(PRE_CONFIG) space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_400_sp_2.csv
	@mv error.csv error_ni_400_sp_2.csv
	@./gridgen 801
	@$(PRE_CONFIG) space_order=1
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_800_sp_1.csv
	@mv error.csv error_ni_800_sp_1.csv
	@$(PRE_CONFIG) space_order=2
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_ni_800_sp_2.csv
	@mv error.csv error_ni_800_sp_2.csv
	@paste -d " " - - < errors.txt > temp.txt
	@mv temp.txt errors.txt
	@gnuplot -p plot_sens.gnuplot
	@gnuplot -p plot_error.gnuplot

rot: compile gridgen str.dx
	@./gridgen 101 0.0
	@$(MAKE) run
	@cp data_1d.csv data_1d_00.csv
	@./gridgen 101 45.0
	@$(MAKE) run
	@cp data_1d.csv data_1d_45.csv
	@./gridgen 101 90.0
	@$(MAKE) run
	@cp data_1d.csv data_1d_90.csv
	@gnuplot -p plot_rot.gnuplot

time: compile gridgen str.dx
	@rm -f errors.txt
	@./gridgen 101
	@$(PRE_CONFIG) iterations=100 sol_out_screen_int=1 timestep=1.0E-03
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_dt_1E-3.csv
	@mv error.csv error_dt_1E-3.csv
	@$(PRE_CONFIG) iterations=1000 sol_out_screen_int=10 timestep=1.0E-04
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_dt_1E-4.csv
	@mv error.csv error_dt_1E-4.csv
	@$(PRE_CONFIG) iterations=10000 sol_out_screen_int=100 timestep=1.0E-05
	@$(PROGRAM)
	@$(EXTRACT_1D)
	@$(CALC_ERROR)
	@mv data_1d.csv data_dt_1E-5.csv
	@mv error.csv error_dt_1E-5.csv
	@gnuplot -p plot_time.gnuplot
	@gnuplot -p plot_error_time.gnuplot

clean:
	@rm -vf *.o *.mod *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.dx stubex *.csv *.h5 *.xdmf fort.*
	@rm -vrf *~
