set terminal pdf
set output "f.pdf"
set title "Evaluation of the unadjusted function for small values."
set logscale x
set grid x,y
set format x "%2.0t{/Symbol \327}10^{%L}"
plot "f.txt" using 1:2 with lines title "f(x)"


set output "adjustedF.pdf"
set title "Evaluation of the adjusted function for small values."
set logscale x
set grid x,y
set yrange [0.35:0.6]
set format x "%2.0t{/Symbol \327}10^{%L}"
plot "adjustedF.txt" using 1:2 with lines title "adjusted f(x)"




