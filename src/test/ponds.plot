! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out
 
set term @PNG enhanced size 1024,512 font ",8"

unset key
set hidden3d
unset xtics
unset ytics
unset border
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
    	              0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
		      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
		      0.875 0.9333 0 0, 1 0.498 0 0 )

dry=1e-3
set view 28,55

set multiplot layout 2,2 scale 1,1.3
splot './eta-1' u 1:2:4 every 3:3 w l, \
      './eta-1' u 1:2:($3>dry?$3+$4:1e1000):(sqrt($5**2+$6**2)) w pm3d
splot './eta-2' u 1:2:4 every 3:3 w l, \
      './eta-2' u 1:2:($3>dry?$3+$4:1e1000):(sqrt($5**2+$6**2)) w pm3d
splot './eta-3' u 1:2:4 every 3:3 w l, \
      './eta-3' u 1:2:($3>dry?$3+$4:1e1000):(sqrt($5**2+$6**2)) w pm3d
splot './eta-8' u 1:2:4 every 3:3 w l, \
      './eta-8' u 1:2:($3>dry?$3+$4:1e1000):(sqrt($5**2+$6**2)) w pm3d
unset multiplot

set output 'level.png'

dry = 0.
set multiplot layout 2,2 scale 1,1.3
splot './level-1' u 1:2:4 every 3:3 w l, \
      './level-1' u 1:2:($3>dry?$3+$4:1e1000):5 w pm3d
splot './level-2' u 1:2:4 every 3:3 w l, \
      './level-2' u 1:2:($3>dry?$3+$4:1e1000):5 w pm3d
splot './level-3' u 1:2:4 every 3:3 w l, \
      './level-3' u 1:2:($3>dry?$3+$4:1e1000):5 w pm3d
splot './level-8' u 1:2:4 every 3:3 w l, \
      './level-8' u 1:2:($3>dry?$3+$4:1e1000):5 w pm3d
unset multiplot

set output 'vectors.png'

reset
set xrange [0:1000]
set yrange [0:1000]
set size ratio -1
unset key
dry=1e-3
plot './eta-8' u 1:2:($3>dry?$5*30.:0):($3>dry?$6*30.:0) every 2:2 w vec lc 0

! rm -f eta-? level-?
