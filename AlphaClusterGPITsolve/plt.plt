unset multiplot; reset

set size ratio 11/6
set key outside

set xlabel 'r[fm]'
set ylabel '��(r)[fm^{-3/2}]'

plot "gp.txt" w l ti "��(r)"
pause -1