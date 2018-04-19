set palette defined ( 0. "#ff1428", 1/6. "#ff2f92", 1/3. "#ff8ad8", 1/2. "#3791e6", 2/3. "#96e650", 5/6. "#fffb00", 1. "#ffffff")
#set palette defined ( 0. "#080818", 1/6. "#012093", 1/3. "#0433ff", 1/2. "#b628b4", 2/3. "#da225c", 5/6. "#ff9300", 1. "#ffcc00")
set size ratio -1
set view map
a=1
b=1
xc=0
yc=0
xmin=xc-a;xmax=xc+a;ymin=yc-b;ymax=yc+b
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 800,700 font 'Helvetica, 14'
set lmargin at screen 0.1
set rmargin at screen 0.86
set tmargin at screen 0.96
set bmargin at screen 0.12
set output 'OUTFILE.png'
HEA:set label 1 't = TIME' left at screen 0.085,0.05 tc rgbcolor "#0000ff" front
set contour
unset surface
set view map
set cntrparam levels discrete 0
#set cntrparam levels incremental -0.8,0.2,0.8
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
set output
