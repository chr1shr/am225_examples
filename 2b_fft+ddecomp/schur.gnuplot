set terminal png size 800,400
set lmargin at screen 0.1
set rmargin at screen 0.85
set output 'schur.png'
set autoscale fix
set pm3d map
splot 'out' matrix using ($1/100.0-1):(1-$2/100.0):3
