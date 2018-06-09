set term pdf color enhanced font "Helvetica"
set auto
set border 15 lw 1

set xlabel "{/Symbol=\326}s [GeV]"

set out './output/cross-section.pdf'
set multiplot
  p "./output/GKPRY.dat" u 1:6 t "Real" with l lw 3 lc "red"
unset multiplot
