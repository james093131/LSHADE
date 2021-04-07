reset
set term gif animate delay 0.05
#set term gif animate
set output "R50.gif"
# set xrange [-0.2:0.2]
set xrange [-5:10]
set yrange [-5:10]
# set yrange [-0.2:0.2]
set zrange [0:1000.0]
set xyplane at 0

# do for [i=0:1999:10]{
do for [i=0:4900:100]{
  set title sprintf("RASTRIGIN time=%i",i)
  #splot sprintf("./A2/A2%i.txt",i) w  lp  t "net"
  splot sprintf("./R50/R50%i.txt",i) w  lp  t "net"
  #splot sprintf("./B30/B30%i.txt",i) w lp t "net", sprintf("test/ball_radviz/ball_radviz_%i.txt",i) w p t "ball"
}
set output
