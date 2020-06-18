set size ratio -1
set key outside
v=0
plot 'cells' w l lc 0, \
  'stencil' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v), \
  'coarse' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''
