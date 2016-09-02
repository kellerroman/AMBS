FC=gfortran

HOMEDIR=../../

CGNSlibs = /home/kell_ro/software/cgns/src/cgns_build/lib/libcgns.a /media/LDS/module/tools/hdf5/1.8.15_GCC-4.9.2/lib/libhdf5.a
LIBS = $(CGNSlibs) -lm -ldl -lz

INCLUDE = -I/usr/include/ -I/home/kell_ro/software/cgns/src/cgns_build/include/
CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries -Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

PRE_CONFIG=$(HOMEDIR)tools/bin/pre_config
#SOL2VTK=$(HOMEDIR)tools/bin/sol2vtk
SOL2VTK=/usr/bin/cgnstools/cgns_to_vtk data_out.cgns
PROGRAM=$(HOMEDIR)bin/AMBS

.PHONY: all cmp run clean cp grid_inf show

	
all: cmp run #show

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

	

run: data_in.cgns config.bin
	@$(PROGRAM) #-debug
	@$(SOL2VTK)
show:
	@gnuplot -p plot.gnuplot

data_in.cgns: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)
gridgen: gridgen.F90
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
