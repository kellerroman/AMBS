FC=gfortran

HOMEDIR=../../

CGNSlibs = -lcgns -lhdf5
LIBS = $(CGNSlibs) -lm

INCLUDE = -I/usr/include/
CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries -Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

PRE_CONFIG=$(HOMEDIR)tools/bin/pre_config
#SOL2VTK=$(HOMEDIR)tools/bin/sol2vtk
SOL2VTK=/usr/bin/cgnstools/cgns_to_vtk data_out.cgns
PROGRAM=$(HOMEDIR)bin/AMBS

.PHONY: all cmp run clean cp grid_inf show

	
all: cmp run post #show

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

	

run: data_in.cgns config.bin
	@$(PROGRAM) #-debug
post: flate_plate_extract
	@$(SOL2VTK)
	@./flate_plate_extract
show:
	@gnuplot -p plot.gnuplot

data_in.cgns: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)

gridgen: gridgen.F90
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

flate_plate_extract:flate_plate_extract.F90 const_mod.o control_mod.o data_mod.o file_io_mod.o screen_io_mod.o
	$(FC) $(CFLAGS) $^ $(LIBS) -o $@ $(INCLUDE)

%.o:../../src/%.F90
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
const_mod.o:../../src/const_mod.F90
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
data_mod.o:../../src/data_mod.F90
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
file_io_mod.o:../../src/file_io_mod.F90 screen_io_mod.o
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
screen_io_mod.o:../../src/screen_io_mod.F90
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
control_mod.o:../../src/control_mod.F90
	$(FC) $(CFLAGS) $< -c $(INCLUDE)
clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
