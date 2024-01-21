set terminal png size 800,600
set title 'NBody3D fifth version performance'
set output '../results/perf_plot/nbody5.png'
set style data histograms
set style fill solid 1.00 border -1
set auto x
set yrange [0:*]
set grid y
set ylabel "GFLOP/s"
set xlabel "Compilation Flags"
plot '../results/simplified_data/gcc/data_5.dat' using 2:xtic(1) title 'GCC' linecolor rgb "blue", \
     '../results/simplified_data/clang/data_5.dat' using 2:xtic(1) title 'CLANG' linecolor rgb "red"
