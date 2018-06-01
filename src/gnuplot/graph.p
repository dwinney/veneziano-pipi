set term pdf color enhanced font "Helvetica"
set auto
set border 15 lw 1

### Angular Dependence
set xlabel "s [GeV^{-2}]"
set title "n_{max} = 6"
set out './output/plots/pwave.pdf'
set ylabel "f^{(1)}_1(s)"
set multiplot
p "./output/rho.dat" u 1:2 t "Re" with l lw 3 lc "red", \
  "./output/rho.dat" u 1:3 t "Im" with l lw 3 lc "blue"
