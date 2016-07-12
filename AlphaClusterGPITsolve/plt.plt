unset multiplot; reset

set size ratio 11/6
set key outside

set xlabel 'r[fm]'
set ylabel 'ƒÌ(r)[fm^{-3/2}]'

plot "gp.txt" w l ti "ƒÌ(r)"
pause -1