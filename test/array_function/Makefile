CFLAGS += -Wconversion -Wall -Wextra -Warray-temporaries -Werror
CFLAGS += -O3
CFLAGS += -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -fbounds-check


all: func_test
	./$^


func_test: func_test.F90
	gfortran $(CFLAGS) $^ -o $@


clean:
	@rm -vf func_test *.mod
