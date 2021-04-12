set xlabel "x"
set ylabel "y"
m = "./plot.dat"
set terminal x11 0
set nokey
set grid
set title "My data from Fortran"
plot m using 1:2 with lines
