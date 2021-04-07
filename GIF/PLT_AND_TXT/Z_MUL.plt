

do for [i=0:8:2]{
  reset
  set term gif animate delay 0.05
  #set term gif animate
  set xrange [-5: 10]
  set yrange [-5:10]
  set zrange [0:30000.0]
  set view 60, 50, 1, 1
  # set view 15, 30, 1, 1
    # set view 60, 75, 1, 1
    # set view 60, 15, 1, 1


  set xyplane at 0
  set output sprintf("../FunctionGIF/Z/Z10_%i.gif",i+1)
  do for [j=0:4900:100]{
    set title sprintf("Zakharov time=%i D%i-D%i",j,i+1,i+2)
    splot sprintf("Z10/Z10%i.txt",j) u i+1:i+2:11 ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  }
  set output
}
  
