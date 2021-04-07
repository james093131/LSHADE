

  reset
  set term gif animate delay 0.05
  #set term gif animate
  set xrange [-5: 10]
  set yrange [-5:10]
  set zrange [0:2000.0]
  # set view 60, 50, 1, 1
  # set view 15, 30, 1, 1
    # set view 60, 15, 1, 1
  set view 65, 330, 1, 1


  set xyplane at 0
  set output sprintf("../FunctionGIF/Z/Z10_DE.gif")
  do for [j=0:5000:100]{
    set title sprintf("Zakharov time=%i ",j)
    splot sprintf("/home/ailab/Downloads/LSHADE/RECORD/Z10_%i.txt",j) u 1:2:11 ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  }
  set output

  
