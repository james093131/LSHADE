do for [i=0:8:2]{
  reset
  set term gif animate delay 0.05
  #set term gif animate
  set xrange [-100:100]
  set yrange [-100:100]
  set zrange [0:30000000000.0]
  set xyplane at 0
  set view 60, 135, 1, 1
  set output sprintf("../FunctionGIF/B/B10_%i.gif",i+1)
  do for [j=0:4900:100]{
    set title sprintf("Bent Cigar time=%i D%i-D%i",j,i+1,i+2)
    splot sprintf("B10/B10%i.txt",j) u i+1:i+2:11 ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  }
  set output
}
  
