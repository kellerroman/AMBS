echo "Parameter Variation for the DIHT Test Case"
HOMEDIR=../../..
PRE=$HOMEDIR/tools/bin/pre_config
SOLVER=$HOMEDIR/bin/AMBS
POST=DIT_gridgen/FFT/fft
make cmp

$PRE disc/space_order=1
$SOLVER
$POST
mv spectrum.dat spectrum_1st_order.dat


$PRE disc/space_order=2
$SOLVER
$POST
mv spectrum.dat spectrum_2nd_order.dat

$PRE disc/space_order=3
$SOLVER
$POST
mv spectrum.dat spectrum_3rd_order.dat


echo 'set logscale x'  > plot_para_var.gnuplot
echo 'set logscale y' >> plot_para_var.gnuplot
echo 'set yrange [1:1000]' >> plot_para_var.gnuplot
echo 'set xrange [0.3:10]' >> plot_para_var.gnuplot
echo 'plot "../data/measurements.dat" u 1:2 \' >> plot_para_var.gnuplot
echo '   , ""                         u 1:3 \' >> plot_para_var.gnuplot
echo '   , ""                         u 1:4 \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:2  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:3  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:4  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:5  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:6  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:7  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:8  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:9  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:10 w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_1st_order.dat"   u 1:11 w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:2  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:3  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:4  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:5  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:6  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:7  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:8  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:9  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:10 w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_2nd_order.dat"   u 1:11 w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:2  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:3  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:4  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:5  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:6  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:7  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:8  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:9  w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:10 w l \' >> plot_para_var.gnuplot
echo '   , "spectrum_3rd_order.dat"   u 1:11 w l ' >> plot_para_var.gnuplot
gnuplot -preserve plot_para_var.gnuplot
