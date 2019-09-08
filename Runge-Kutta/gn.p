set xlabel 'dt'
set ylabel 'error'
set logscale
set key top left
set ytics format '%.0e'
plot "< grep '1$' log" pt 7 t '', 15.*x t '15 dt',  \
     "< grep '2$' log" pt 7 t '', 4.*x*x t '4 dt^2', \
     "< grep '4$' log" pt 7 t '', x**4/2. t 'dt^4/2'
