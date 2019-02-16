plot [-0.5:0.5]'../yprof.ghia' u 1:2 title "Ghia et al." w p ps 1.5 pt 6 lw 2, \
     'yprof' w l lw 2 title "Basilisk (centered)", \
     '../lidmac/yprof' w l lw 2 title "Basilisk (MAC)"
