FC=gfortran

HOMEDIR=../../
CODE_OBJ=$(HOMEDIR)obj/
TOOLS_OBJ=$(HOMEDIR)tools/obj/

LIBS = -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lhdf5_fortran -lm

INCLUDE = -I/usr/include/ -I$(TOOLS_OBJ) -I$(CODE_OBJ)
CFLAGS = -std=f2008
CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries -Werror
#CFLAGS += -DDEBUG #Unnecessary comments are oppressed
CFLAGS += -O3

PRE_CONFIG=$(HOMEDIR)tools/bin/pre_config
CREATE_XDMF=$(HOMEDIR)tools/bin/create_xdmf
PROGRAM=$(HOMEDIR)bin/AMBS

.PHONY: all cmp run clean cp grid_inf show

	
all: compile run post #show

compile: compile_source flate_plate_extract
compile_source:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

run: data_in.h5 config.bin
	@$(PROGRAM) #-debug
	@$(CREATE_XDMF)
post: compile
	@./flate_plate_extract
show:
	@gnuplot -p plot.gnuplot

data_in.h5: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

flate_plate_extract: flate_plate_extract.F90 $(TOOLS_OBJ)module_post.o $(CODE_OBJ)const_mod.o $(CODE_OBJ)file_io_mod.o $(CODE_OBJ)screen_io_mod.o $(CODE_OBJ)data_mod.o $(CODE_OBJ)control_mod.o
	$(FC) $(CFLAGS) $^ $(LIBS) -o $@ $(INCLUDE)

clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.dx *.h5 *.csv flate_plate_extract *.xdmf fort.*
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
