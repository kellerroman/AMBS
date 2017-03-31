FC=gfortran

HOMEDIR=../../../
CODE_OBJ=$(HOMEDIR)obj/
TOOLS_OBJ=$(HOMEDIR)tools/obj/

LIBS = -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lhdf5_fortran -lm -lfftw3

INCLUDE = -I/usr/include/ -I$(TOOLS_OBJ) -I$(CODE_OBJ)
#CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries #-Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

PRE_CONFIG=$(HOMEDIR)tools/bin/pre_config
CREATE_XDMF=$(HOMEDIR)tools/bin/create_xdmf
PROGRAM=$(HOMEDIR)bin/AMBS

.PHONY: all compile run clean cp show

	
all: compile run post show
#all: cmp export_velocities run post exp_vel show

exp_vel: export_velocities
	@./$<
	@DIT_gridgen/FFT/fft

compile: compile_source export_velocities

compile_source:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory


run: compile data_in.h5 config.bin
	@$(PROGRAM)
post: 
	@$(CREATE_XDMF)
show:
	@gnuplot -p plot.gnuplot

data_in.h5: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)
	
gridgen: gridgen.F90 $(TOOLS_OBJ)module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)
	
export_velocities: export_velocities.F90 $(TOOLS_OBJ)module_post.a 
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

clean:
	@rm -vf bin *.o *.mod *.h5 *.dat *.plt *.png *.xdmf *.txt *.bin gridgen export_velocities *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
