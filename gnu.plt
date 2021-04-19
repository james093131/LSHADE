set terminal png
set output "DIM_30_A_PART.png"
set xrange [1000:10000]

plot "ACKLEY/DIM_30_A4.txt" w lp ps 0.8  lw 3 lc rgb '#778da9' t "line"

set output
