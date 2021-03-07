reset
set title 'Time-reversed advection'

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set xrange [32:512]
set xtics 32,2,512
set grid ytics
set cbrange [1:2]
plot 'log' u 1:4 t 'max', 'log' u 1:2 t 'norm1', exp(f(log(x))) t ftitle(a,b), exp(f2(log(x))) t ftitle(a2,b2)

if (batch) set term @PNG; set output "error.png"; else pause -1;
reset
set title 'Time-reversed advection'
set pm3d
set pm3d map interpolate 1,1
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,\
     0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,\
     0.875 0.9333 0 0, 1 0.498 0 0 )
set size ratio -1
splot 'out' t ""
