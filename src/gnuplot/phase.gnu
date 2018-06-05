set term pdf color enhanced font "Helvetica"
set auto
set border 15 lw 1

set xlabel "{/Symbol=\326}s [GeV]"

set ylabel "{/Symbol d}(s)"
set out './output/phase.pdf'
  p "./output/phase_shift.dat" u 1:2 t "S0" with l lw 3 lc "red"

set ylabel "{/Symbol h}(s)"
set out './output/inelasticty.pdf'
  p "./output/phase_shift.dat" u 1:3 t "S0" with l lw 3 lc "red"
