#set datafile separator ','
#set grid
# set term postscript eps enhanced color
#set term png
set ylabel "Y (m)"
set xlabel "Geschwindigkeit"
#set logscale y
set title "Profile"
#set output 'rms.png'
plot 'profile.csv' using 2:1 notitle with lines\
    , ""  u 3:1 notitle with lines\
    , ""  u 4:1 notitle with lines\
    , ""  u 5:1 notitle with lines\
    , ""  u 6:1 notitle with lines\
    , ""  u 7:1 notitle with lines\
    , ""  u 8:1 notitle with lines\
    , ""  u 9:1 notitle with lines\
    , ""  u 10:1 notitle with lines\
    , ""  u 11:1 notitle with lines\
    , ""  u 12:1 notitle with lines\
    , ""  u 13:1 notitle with lines\
    , ""  u 14:1 notitle with lines\
    , ""  u 20:1 \
    , ""  u 30:1 notitle with lines\
    , ""  u 40:1 notitle with lines\
    , ""  u 50:1 notitle with lines\
    , ""  u 60:1 notitle with lines\
    , ""  u 100:1 notitle with lines
#, 'RMS.csv' using 1:2 notitle with lines
