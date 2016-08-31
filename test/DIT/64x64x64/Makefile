FC=gfortran

HOMEDIR=../../../

CGNSlibs = -lcgns -lhdf5
LIBS = $(CGNSlibs) -lm

INCLUDE = -I/usr/include/
CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries -Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

PRE_CONFIG=$(HOMEDIR)tools/bin/pre_config
EXTRACT_1D=$(HOMEDIR)tools/bin/extract_1D
SOL2VTK=/usr/bin/cgnstools/cgns_to_vtk data_out.cgns
PROGRAM=$(HOMEDIR)bin/AMBS

.PHONY: all cmp run clean cp grid_inf show

	
all: cmp run post exp_vel

exp_vel: export_velocities
	@./$<
	@DIT_gridgen/FFT/fft
cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory


run: data_in.cgns config.bin
	@$(PROGRAM)
post: 
	@$(SOL2VTK)
show:
	@gnuplot -p plot.gnuplot

data_in.cgns: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)
	
gridgen: gridgen.F90 
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)
	
export_velocities: export_velocities.F90 file_io_mod.o const_mod.o data_mod.o screen_io_mod.o control_mod.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

const_mod.o:../../../src/const_mod.F90
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
data_mod.o:../../../src/data_mod.F90
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
file_io_mod.o:../../../src/file_io_mod.F90 screen_io_mod.o data_mod.o
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
screen_io_mod.o:../../../src/screen_io_mod.F90 control_mod.o
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
control_mod.o:../../../src/control_mod.F90
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
