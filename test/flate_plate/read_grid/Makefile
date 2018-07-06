HOMEDIR=../../../

include $(HOMEDIR)Makefile.global

INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

.PHONY: all cmp run clean cp grid_inf show
.DEFAULT_GOAL := all

	
all: compile run post #show

compile: compile_global flate_plate_extract

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

gridgen: gridgen.F90 $(TOOLS_OBJ)/module_gridgen.o $(CODE_OBJ)const_mod.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)

flate_plate_extract: flate_plate_extract.F90 $(TOOLS_OBJ)module_post.o $(CODE_OBJ)const_mod.o $(CODE_OBJ)file_io_mod.o $(CODE_OBJ)screen_io_mod.o $(CODE_OBJ)data_mod.o $(CODE_OBJ)control_mod.o
	$(FC) $(CFLAGS) $^ $(LIBS) -o $@ $(INCLUDE)

clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.dx *.h5 *.csv flate_plate_extract *.xdmf fort.*
	@rm -vrf *~
