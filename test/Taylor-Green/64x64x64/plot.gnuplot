set logscale x
set logscale y
set yrange [1:1000]
set xrange [0.3:10]
plot "../data/measurements.dat" u 1:2 \
   , ""                         u 1:3 \
   , ""                         u 1:4 \
   , "spectrum_1000.dat" u 1:2  w l \
   , "spectrum.dat"      u 1:2  w l \
   , "spectrum.dat"      u 1:3  w l \
   , "spectrum.dat"      u 1:4  w l \
   , "spectrum.dat"      u 1:5  w l \
   , "spectrum.dat"      u 1:6  w l \
   , "spectrum.dat"      u 1:7  w l \
   , "spectrum.dat"      u 1:8  w l \
   , "spectrum.dat"      u 1:9  w l \
   , "spectrum.dat"      u 1:11 w l \
   , "spectrum.dat"      u 1:11 w l 

