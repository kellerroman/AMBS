HOMEDIR=../../../

include $(HOMEDIR)/Makefile.global
INCLUDE += -I$(TOOLS_OBJ) -I$(CODE_OBJ)

.PHONY: all run clean show
.DEFAULT_GOAL := all

	
all: compile run post show

exp_vel: export_velocities
	@./$<

compile: compile_global export_velocities

run: compile data_in.h5 config.bin
	@$(PROGRAM)

post: export_velocities 
	@$(CREATE_XDMF)
	@./$<
show:
	#@gnuplot -p plot.gnuplot

data_in.h5: gridgen
	@./$<

config.bin: config.cfg	
	@$(PRE_CONFIG)
	
gridgen: gridgen.F90 $(TOOLS_OBJ)module_gridgen.o
	$(FC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)
	
export_velocities: export_velocities.F90 $(TOOLS_OBJ)module_post.a 
	#$(FC) $(CFLAGS) -Wno-unused-parameter  $^ -o $@ $(LIBS) $(INCLUDE)

clean:
	@rm -vf bin *.o *.mod *.h5 *.dat *.plt *.png *.xdmf *.txt *.bin gridgen export_velocities *.silo *.dx stubex *.cgns *.csv
	@rm -vrf *~
