

  reset
  set term gif animate delay 0.05
  #set term gif animate
  set xrange [-100: 100]
  set yrange [-100:100]
  set zrange [0.0:5.0]
  set view 40, 30, 1, 1
  # set view 15, 30, 1, 1
    # set view 60, 75, 1, 1
    set view 65, 35, 1, 1


  set xyplane at 0
  set output sprintf("../FunctionGIF/DE/S10_DE.gif")
  do for [j=0:5000:100]{
    set title sprintf("Schaffer_F7 time=%i ",j)
    splot sprintf("/home/ailab/Downloads/LSHADE/RECORD/S10_%i.txt",j) u 1:2:11 ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  }
  set output

  