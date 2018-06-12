set term pdf color enhanced font "Helvetica"
set auto
set border 15 lw 1

set xlabel "{/Symbol=\326}s [GeV]"
set xrange[.0:2.1]
set out './output/fit.pdf'

set multiplot
p "./output/fit.dat" u 1:2 t "Venez. Real" with l lw 3 lc "red", \
  "./output/fit.dat" u 1:3 t "Venez. Imaginary" with l lw 3 lc "blue",\
  "./output/data.dat" u 1:2 t "GKPY Real" w l lc "red" dt 3 lw 3, \
  "./output/data.dat" u 1:3 t "GKPY Imaginary" w l lc "blue" dt 3 lw 3
