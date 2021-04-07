reset
set term gif animate delay 0.05
#set term gif animate

set output "../FunctionGIF/DE/R10_DE.gif"
set xrange [-5.12: 5.12]
set yrange [-5.12:5.12]
set zrange [0.0:150.0]
set xyplane at 0

do for [i=0:5000:100]{
  set title sprintf("RASTRIGIN time=%i",i)
  splot sprintf("/home/ailab/Downloads/LSHADE/RECORD/R10_%i.txt",i) u 1:2:11 ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  #bfa89e
}
set output
