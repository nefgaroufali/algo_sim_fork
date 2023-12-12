# plot_script2.gp

set term x11 persist # Use X11 terminal for displaying the plot (you can choose other terminals as well)
set title "Your Plot Title"
set xlabel "X-axis label"
set ylabel "Y-axis label"

plot 'file_1.txt' using 1:2 with linespoints title "Data Series 1"