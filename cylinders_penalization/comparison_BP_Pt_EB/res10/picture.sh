#!/usr/bin/gnuplot -c
# Show usage information if no argument is present
if (strlen(ARG1) == 0) print "Usage: " . ARG0 . " level"; exit
# Output W3C Scalable Vector Graphics
set terminal postscript eps enhanced color font 'Times-Roman,10'
level=ARG1
set output "comparison10.eps"
array(n) = word("Pt10 1e-15 EB10", n)
titles(n) = word("'Simple Splitting method (SSM)' 'Brinkman Penalization method (BPM), eta=1e-15' 'Embedded Boundary method (EBM)' 'Theory: Sangani and Acrivos, 1982'",n)
N=3
grepf(l, num) = sprintf("< grep \"^%s\" log_%s", l, num)
set bmargin 5
set lmargin 13
set xlabel 'Volume fraction'  font ',18'  offset 0,-1,0
set ylabel 'k'  font ',18' offset -4,0,0
set key font ",15"
set logscale y
set grid
set key spacing 2
set key top left
set tics font ",16"
set format y "10^{%L}"
set ytics 10
set xtics 0.1
array files[N+1]
array lstyles[N+1]
array col[N+1]
set for [i=1:7] linetype i dt i
#set style line 1 lt 1 lc rgb "blue" lw 3 pt 6 ps 1
#set style line 2 lt 1 lc rgb "orange" lw 3 pt 6 ps 1
set style line 1 lt 1 lc rgb "gold" lw 3 pt 30 ps 1
#set style line 4 lt 3 lc rgb "green" lw 3 pt 33 ps 1
set style line 2 lt 1 lc rgb "blue" lw 3 pt 10 ps 1
set style line 3 lt 2 lc rgb "red" lw 3 pt 12 ps 1
set style line 4 lt 3 lc rgb "black" lw 3 pt 2 ps 1

#show style line
p=1
do for [i=1:N+1] {
	files[p]=grepf(level, array(i))
	lstyles[p]=i
	col[p]=3
	p=p+1
}
files[N+1]=grepf(level, array(1))
col[N+1]=4
plot for [i=1:N+1] files[i] u 2:col[i] w lp ls lstyles[i] t titles(i)
