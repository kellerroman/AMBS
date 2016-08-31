EXECUTABLE=bfs_gridgen

#CFLAGS+= -Wall                     # Enable most warning messages.
                                    # Enables commonly used warning options pertaining to 
                                    # usage that we recommend avoiding and that we believe 
                                    # are easy to avoid. 
#CFLAGS+= -Wextra                   # Enable warning not enabled by -Wall.
#CFLAGS+= -Wuninitialized           # Warn about variables used before initialized. 
                                    # Requires -O2. Only limited practical use.
                                    # Disabled if -finit is activated.
 CFLAGS+= -fimplicit-none           # Do not allow implicit declaration.
 CFLAGS+= -fbounds-check            # Check array and substring bounds at runtime. 
 CFLAGS+= -fbacktrace               # Will insert code to allow stack tracebacks on abnormal 
                                    # end of program. This will slow down your program.
 CFLAGS+= -finit-integer=-12345     # Initialize uninitialized scalar integer variables to n: 
                                    # -finit-integer=n.
 CFLAGS+= -finit-real=nan           # Initialize uninitialized scalar real variables.
                                    # Legal values are none,zero, nan, inf, +inf and -inf.
                                    # Finally, note that enabling any of the ‘-finit-*’ i
                                    # options will silence  warnings that would have been 
                                    # emitted by ‘-Wuninitialized’ for the  affected local 
                                    # variables.
 CFLAGS+= -ffpe-trap=invalid,zero,overflow # stop at first occurence of invalid number
 CFLAGS+= -g
all: git.bin


git.bin: $(EXECUTABLE)
	./$<

$(EXECUTABLE):$(EXECUTABLE).F90
	gfortran $(CFLAGS) $< -o $@
plot: git.bin
	u2p -git -s
clean:
	rm -rfv $(EXECUTABLE) git.bin zs.dat git.plt
