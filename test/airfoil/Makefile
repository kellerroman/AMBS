FC=gfortran

HOMEDIR=../../

LIBS = -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lhdf5_fortran -lm

INCLUDE = -I/usr/include/ -I$(HOMEDIR)tools/obj/ -I$(HOMEDIR)obj/
CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries -Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

PRE_CONFIG=$(HOMEDIR)tools/bin/pre_config
PROGRAM=$(HOMEDIR)bin/AMBS

.PHONY: all cmp run clean cp grid_inf show

	
all: cmp run #show

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

	

#run: gridgen
#	@./gridgen

run: data_in.h5 config.bin
	@$(PROGRAM) #-debug

show:
	@gnuplot -p plot.gnuplot

data_in.h5: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)
gridgen: gridgen.F90 $(HOMEDIR)tools/obj/module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

clean:
	@rm -vf *.o *.mod *.ufo *.dat *.plt *.png *.txt *.bin gridgen *.silo *.dx *.csv *.h5 *.xdmf
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
