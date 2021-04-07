


  reset
  set term gif animate delay 0.05
  #set term gif animate
  set xrange [-5:10]
  set yrange [-5:10]
  set zrange [0:100000.0]
  # set view 50, 325, 1, 1
  set xyplane at 0
  set output sprintf("../FunctionGIF/DE/RO10_DE.gif")
  do for [j=0:5000:100]{
    set title sprintf("Rosenbrock time=%i",j)
    splot sprintf("/home/ailab/Downloads/LSHADE/RECORD/RO10_%i.txt",j) u 1:2:11  ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  }
  set output

