set term pdf color enhanced font "Helvetica"
set auto
set border 15 lw 1

set xlabel "{/Symbol=\326}s [GeV]"

set out './output/rho.pdf'

set multiplot
p "./output/rho.dat" u 1:2 t "S0" with l lw 3 lc "red", \
#  "./output/rho.dat" u 1:3 t "Im" with l lw 3 lc "blue"