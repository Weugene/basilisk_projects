#!/usr/bin/env gnuplot -c
# Show usage information if no argument is present
if (strlen(ARG1) == 0) print "Usage: " . ARG0 . " level"; exit
# Output W3C Scalable Vector Graphics
set terminal postscript eps enhanced color font 'Helvetica,10'
level=ARG1
emin=2
emax=6
N=emax-emin+1
grepf(l, num) = sprintf("< grep \"^%s\" log_1e-%d", l, num)
set xlabel 'Volume fraction'  font ",10"
set ylabel 'k_0'  font ",10"
set logscale y
set grid
set key spacing 2
set key top left
set tics font "Helvetica,10"
set format y "10^{%L}"
set ytics 100
set xtics 1
array files[N+1]
array titles[N+1]
array lstyles[N+1]
set for [i=1:7] linetype i dt i
set style line 1 lt 1 lc rgb "blue" lw 3 pt 6 ps 1
set style line 2 lt 1 lc rgb "orange" lw 3 pt 6 ps 1
set style line 3 lt 1 lc rgb "gold" lw 3 pt 30 ps 1
set style line 4 lt 1 lc rgb "green" lw 3 pt 33 ps 1
set style line 5 lt 1 lc rgb "blue" lw 3 pt 10 ps 1
set style line 6 lt 1 lc rgb "red" lw 3 pt 12 ps 1
set style line 7 lt 3 lc rgb "black" lw 3 pt 2 ps 1

#show style line
p=1
do for [i=emin:emax+1] {
	files[p]=grepf(level, i)
	titles[p]=sprintf("BP eta=1e-%d, %s levels", i, level)
	lstyles[p]=i
	p=p+1
}
files[N+1]=grepf(level, emax)
titles[N+1]=sprintf("Sangani and Acrivos, 1982")
plot for [i=1:N+1] files[i] u 2:3 w lp ls lstyles[i] t titles[i]
