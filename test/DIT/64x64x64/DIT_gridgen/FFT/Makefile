all:
	gcc -o fft fft.c -I./FFTW/fftw2_libs/include/ -I/usr/include -L ./FFTW/fftw2_libs/lib/ -L /usr/lib/ -lrfftw -lfftw -lm -lcgns -lhdf5
clean:
	rm fft
