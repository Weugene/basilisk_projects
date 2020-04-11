#!/usr/bin/gnuplot -c
print "args:".ARG1
if (strlen(ARG1) == 0) print "Usage: " . ARG0 . " picture case. 1-relative Err. 2-average velocity, 3-max velocity"; exit
# Output W3C Scalable Vector Graphics
#set terminal postscript eps enhanced color font 'Helvetica,10'
#set terminal pdf
set terminal postscript eps enhanced color font 'Helvetica,10'
set grid
set key spacing 2
set key top left
set tics font "Helvetica,10"
set for [i=1:7] linetype i dt i
set style line 1 lt 1 lc rgb "blue" lw 3 pt 1 ps 2
set style line 2 lt 1 lc rgb "red" lw 3 pt 2 ps 2
set style line 3 lt 1 lc rgb "green" lw 3 pt 4 ps 2
set style line 4 lt 1 lc rgb "green" lw 3 pt 33 ps 1
set style line 5 lt 1 lc rgb "blue" lw 3 pt 10 ps 1
set style line 6 lt 1 lc rgb "red" lw 3 pt 12 ps 1
set style line 7 lt 3 lc rgb "black" lw 3 pt 2 ps 1
mcase=ARG1+0
array A[3]
array N[3]
array Title[3]
array MAXFi[3]
array lstyles[3]
A[1]=0.08
A[2]=0.125
A[3]=0.2
tmax=0.4
dataname(n) = sprintf("m%g",n)
set xlabel 'time'  font ",10"

do for [i=1:3] {
    N[i]=sprintf("m%g",A[i]);
    Title[i]=sprintf("r1=%g",A[i]);
    stats  N[i] u 3:5 name "XX";
    MAXFi[i]=XX_max_y;
    lstyles[i]=i;
}
set xr [0:tmax]

if (mcase==1){
    set yr [-1e-5:1e-5];
    set format y "10^{%L}"
    set ylabel 'Relative Error, (fi-fi0)/fi0'  font ",10"
    plot for[i=1:3] N[i] u 3:(($5 - MAXFi[i])/MAXFi[i]) t Title[i] w lp ls i;
}

if (mcase==2){
    set ylabel 'Average velocity'  font ",10"
    plot for[i=1:3] N[i] u 3:6 t Title[i]  w lp ls i
}

if (mcase==3){
    set ylabel 'Max velocity'  font ",10"
    plot for[i=1:3] N[i] u 3:8 t Title[i] w lp ls i
}
