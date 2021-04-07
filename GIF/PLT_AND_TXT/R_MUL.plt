do for [i=0:8:2]{
  reset
  set term gif animate delay 0.05
  #set term gif animate
  set xrange [-5.12: 5.12]
  set yrange [-5.12:5.12]
  set zrange [0.0:200.0]
  set xyplane at 0
  set output sprintf("../FunctionGIF/R/R10_%i.gif",i+1)
  do for [j=0:4900:100]{
    set title sprintf("RASTRIGIN time=%i D%i-D%i",j,i+1,i+2)
    splot sprintf("R10/R10%i.txt",j) u i+1:i+2:11 ps 0.2 pt 14  lc rgb '#778da9' t "net" 
  }
  set output
}
