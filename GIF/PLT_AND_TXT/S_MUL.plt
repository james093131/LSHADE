

do for [i=0:8:2]{
  reset
  set term gif animate delay 0.05
  #set term gif animate
  set xrange [-100: 100]
  set yrange [-100:100]
  set zrange [0.0:100.0]
  # set view 60, 45, 1, 1
  # set view 15, 30, 1, 1
    # set view 60, 75, 1, 1
    # set view 60, 15, 1, 1


  set xyplane at 0
  set output sprintf("../FunctionGIF/S/S30_%i.gif",i+1)
  do for [j=0:4900:100]{
    set title sprintf("Schaffer_F7 time=%i D%i-D%i",j,i+1,i+2)
    splot sprintf("S30/S30%i.txt",j) u i+1:i+2:11 ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  }
  set output
}
  