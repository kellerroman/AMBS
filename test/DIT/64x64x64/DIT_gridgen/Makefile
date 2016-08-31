EXECUTABLE=grid_cube

all: git.bin


git.bin: $(EXECUTABLE)
	./$<

$(EXECUTABLE):$(EXECUTABLE).F90
	gfortran $< -o $@
plot: git.bin
	u2p -git -s
clean:
	rm -rfv $(EXECUTABLE) git.bin zs.dat git.plt
