set term pdf color enhanced font "Helvetica"
set auto
set border 15 lw 1

set xlabel "{/Symbol=\326}s [GeV]"
set xrange[.2:1.4]
set out './output/rho.pdf'

set multiplot
p "./output/rho.dat" u 1:2 t "Real" with l lw 3 lc "red", \
  "./output/rho.dat" u 1:3 t "Imaginary" with l lw 3 lc "blue",\
   "./output/rho.dat" u 1:4 t "Real" with l lw 3 lc "magenta", \
    "./output/rho.dat" u 1:5 t "Imaginary" with l lw 3 lc "cyan"
