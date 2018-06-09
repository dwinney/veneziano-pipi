set term pdf color enhanced font "Helvetica"
set auto
set border 15 lw 1

set xlabel "{/Symbol=\326}s [GeV]"

set ylabel "{/Symbol d}(s)"
set out './output/phaseshift.pdf'
  p "./output/GKPRY.dat" u 1:2 t "S0" with l lw 3 lc "red"

set ylabel "{/Symbol h}(s)"
set out './output/inelasticity.pdf'
  p "./output/GKPRY.dat" u 1:3 t "S0" with l lw 3 lc "red"

unset ylabel
set out './output/amplitude.pdf'
set multiplot
  p "./output/GKPRY.dat" u 1:4 t "Real" with l lw 3 lc "red", \
    "./output/GKPRY.dat" u 1:5 t "Imaginary" with l lw 3 lc "blue"
unset multiplot
