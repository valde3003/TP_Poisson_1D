# Set the output format to png
set terminal png size 800,600

# Set the output file name
set output 'Alpha_richardson.png'

# Set titles and labels
set title "RESVEC_RICHARDSON"
set xlabel "Nb éléments"
set ylabel "Temps"

# Assuming a constant time step for each data point
# User can modify this time step as needed
set xtics 1
set grid

# Plotting the data
plot 'RESVEC_RICHARDSON.dat' using ($0):1 with lines title 'ALPHA RICHARDSON'