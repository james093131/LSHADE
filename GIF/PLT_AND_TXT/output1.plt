reset
set term gif animate delay 0.05
#set term gif animate

set output "../FunctionGIF/R2.gif"
set xrange [-5.12: 5.12]
set yrange [-5.12:5.12]
set zrange [0.0:50.0]
set xyplane at 0


  set title "RASTRIGIN "
  splot "R2/R20.txt"  w  lp  t "net"

set output
