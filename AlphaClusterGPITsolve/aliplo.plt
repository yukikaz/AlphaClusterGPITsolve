unset multiplot; reset

set size ratio 11/6
set key outside

set xlabel 'x'
set ylabel '|É’|^2'

plot "ali.txt" w l ti "ali"
pause -1