set macro  # Enable macro definition

# Some macro definitions

# Custom colors: hex #RRGGBB
label_color = "#000000"
tic_color = "#000000"
title_color = "#000000"
myblue = "#2196f3"
myred = "#f44336"
mygreen = "#2e7d32"
mypurple = "#9c27b0"
myviolet = "#673ab7"
mybrown = "#795548"
myorange = "#ff9800"

# Width and height of postscript figure in inches
width = 12.5
height = 12.5

# x-axis resolution
set samples 1000

# Line styles.  
# For lines: plot x with lines
# For points: plot x with points

set style line 1 lc rgb 'black' pt 5   # fill square
set style line 2 lc rgb myred pt 7   # circle
set style line 3 lc rgb 'blue' pt 9   # triangle
set style line 4 lc rgb mygreen pt 7   # circle
set style line 5 lc rgb mypurple pt 7   # circle
set style line 6 lc rgb myviolet  pt 7   # circle
set style line 7 lc rgb mybrown  pt 7   # circle
set style line 8 lc rgb myorange  pt 7   # circle
set style line 9 lc rgb myblue pt 7   # circle
set style line 10 lc rgb 'green' pt 6   # open circle
set style line 11 lc rgb 'red' pt 6   # open circle
set style line 12 lc rgb 'gray' pt 4   # open circle
set style line 13 lc rgb 'gold' pt 4   # open circle
set style line 14 lc rgb 'dark-grey' pt 4   # open circle
set style line 15 lc rgb 'dark-cyan' pt 4   # open circle
set style line 16 lc rgb 'purple' pt 4   # open circle
set style line 17 lc rgb 'orchid' pt 4   # open circle
set style line 18 lc rgb 'brown' pt 4   # open circle
set style line 19 lc rgb 'dark-green' pt 4   # open circle
set style line 20 lc rgb 'magenta' pt 4   # open circle

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log time (s)' textcolor rgb tic_color #font "Arial,22"
set ylabel 'Log X' textcolor rgb tic_color #font "Arial,22"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.0    # Size of the plotted points

#set key right top outside    # Place legend inside top 
unset key            # Don't show legend in screen plot (will show in eps)

set timestamp       # Date/time

ds="C++ Asy 150 with PF"
ds = ds.": viktorExtendedProfileSmooth.inp"
set title ds textcolor rgb title_color


# -------- Axis ranges and ticmarks -----------

xlow = -16
xup = -2 
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -14
yup = 0
ytics = 2    # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

#set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------


modsize = 20          # Number independent linestyles defined above
numberCurves = 15   # Number species to be plotted
widdy = 1.0           # Curve linewidths (approx twice linewidth in pts)
dasher = 1            # Dash style for curves (0,1,2,3, ...)

file1 = "plot1.data"  # Data file with mass fractions X

# Loop to plot X for numberCurves isotopes output from 
# explicitMatrix.cpp -> gnu_out/plot1.data.  There are modsize
# line styles defined above. Use the mod operator % to cycle
# through the set of line styles once for every modsize of the 
# numberCurves isotopes to be plotted.

plot for[i=8 : numberCurves+8] file1 using 1:i with lines \
ls ((i-8)%modsize+1) lw widdy dashtype dasher title "(".(i-8).")"

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,20"
set key top right outside font "Arial,15"
set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,22"
set ylabel 'Log X' textcolor rgb tic_color font "Arial,22"

# Plot to postscript file

#unset key
set out "gnuplot_X-loop.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 18
replot               # Plot to postscript file

# Plot to PNG file (un-comment to enable)

#set out "gnuplot_X-loop.png"
## Assume 72 pixels/inch and make bitmap twice as large for display resolution
#set terminal png transparent size 2*width*72, 2*height*72 lw 2
#replot

quit