HOMEDIR=../../
include $(HOMEDIR)Makefile.global

INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)
.PHONY: all run clean grid_inf show
.DEFAULT_GOAL := all

	
all: compile_global run #show

run: data_in.h5 config.bin
	@$(PROGRAM) #-debug

show:
	@gnuplot -p plot.gnuplot

data_in.h5: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)

gridgen: gridgen.F90 $(TOOLS_OBJ)module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

clean:
	@rm -vf *.o *.mod *.ufo *.dat *.plt *.png *.txt *.bin gridgen *.silo *.dx *.csv *.h5 *.xdmf
	@rm -vrf *~

