unset multiplot; reset

set size ratio 11/6
set key outside

set xlabel 'x'
set ylabel '|É’|^2'

plot "hR.txt" w l ti "hR"
pause -1