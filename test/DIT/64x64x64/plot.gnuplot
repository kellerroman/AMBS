set logscale x
set logscale y
plot "../data/measurements.dat" u 1:2 \
   , ""                         u 1:3 \
   , ""                         u 1:4 \
   , "spectrum_1000.dat" u 1:2 \
   , "spectrum.dat"           u 1:2 

