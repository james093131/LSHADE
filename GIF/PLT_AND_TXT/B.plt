
  reset
  set term gif animate delay 0.05
  #set term gif animate
  set xrange [-100:100]
  set yrange [-100:100]
  set zrange [0:10000000000.0]
  set xyplane at 0
  # set view 60, 135, 1, 1
  set output sprintf("../FunctionGIF/DE/B10_DE.gif")
  set view 60, 130, 1, 1

  do for [j=0:5000:100]{
    set title sprintf("Bent Cigar time=%i ",j)
    splot sprintf("/home/ailab/Downloads/LSHADE/RECORD/B10_%i.txt",j) u 1:2:11 ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  }
  set output

  
