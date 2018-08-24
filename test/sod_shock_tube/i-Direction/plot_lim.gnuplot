plot \
   "str.dx"               u 1:2 w l t "Analytical Solution"\
  ,"data_lim=1.csv" u 1:2 w l t "minmod"\
  ,"data_lim=2.csv" u 1:2 w l t "vanLeer"\
  ,"data_lim=3.csv" u 1:2 w l t "superbee"
