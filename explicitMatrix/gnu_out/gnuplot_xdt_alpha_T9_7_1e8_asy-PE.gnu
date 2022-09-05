set macro  # Enable macro definition

# Some macro definitions

# Colors: hex #RRGGBB
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
width = 8.0  # 8.5
height = 5.0   # 4.5

# x-axis resolution
set samples 1000

# Line styles.  
# For lines: plot x with lines ls 1
# For points: plot x with points ls 1

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
set style line 11 lc rgb 'red' pt 4   # open circle
set style line 12 lc rgb 'gray' pt 4   # open circle
set style line 13 lc rgb 'gold' pt 4   # open circle
set style line 14 lc rgb 'orange' pt 4   # open circle

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log t (s)' textcolor rgb tic_color #font "Arial,32"
set ylabel 'Log dt (s)' textcolor rgb tic_color #font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key outside     # Place legend outside top
#unset key            # Don't show legend

set timestamp       # Date/time

ds="Asy-PE and Asy timestepping"
ds = ds.": alpha T9=7 rho=1e8, with PF"
set title ds textcolor rgb title_color



# -------- Axis ranges and ticmarks -----------

xlow = -16
xup = 0
xtics = 2    # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -17 
yup = -1
ytics = 2      # Space between major y ticmarks
minytics = 5   # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

#set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------


# Set data files

file1 = "plot2.data"  # Data: C++ Asy+PE
file2 = "gnufile_alpha_T9_7_1e8_asy_Java_PF.data"  # Java Asy (1e-9)
file3 = "gnufile_alpha_T9_7_1e8_asy_C++_PF.data"  # C++ Asy (1e-8)

# Plot data

plot file1 using 1:5 with lines ls 1 lw 1.5 dashtype 2 title "log 2/Rmax"
replot file3 using 1:2 with lines ls 3 lw 1.5 dashtype 1 title "log dt (asy C++ 1e-8)"
replot file1 using 1:2 with lines ls 2 lw 1.5 dashtype 1 title "log dt (asy+PE C++)"
replot file2 using 1:2 with lines ls 10 lw 1.5 dashtype 1 title "log dt (asy java 1e-9)"

# Plot dt vs t reference lines

replot file1 using 1:( log10((10**$1)*0.1) ) with lines ls 1 lw 0.6 dashtype 2 title "dt=0.1 t"
replot file1 using 1:( log10((10**$1)*0.01) ) with lines ls 1 lw 0.6 dashtype 0 title "dt=0.01 t"
replot file1 using 1:( log10((10**$1)*0.001) ) with lines ls 1 lw 0.4 dashtype 7 title "dt=0.001 t"
replot file1 using 1:( log10((10**$1)*0.0001) ) with lines ls 1 lw 0.4 dashtype 9 title "dt=0.0001 t"

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,18"
set key top left inside font "Arial,16"
set xlabel 'Log t (s)' textcolor rgb tic_color font "Arial,24"
set ylabel 'Log dt (s)' textcolor rgb tic_color font "Arial,24"

# Plot to postscript file. Use Symbol font because minus signs are 
# lost for many other fonts

set out "gnuplot_xdt_alpha_T9_7_1e8_asy-PE.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Symbol" 22
replot               # Plot to postscript file

# Plot to PNG file

#set out "gnuplot_xdt_alpha_T9_7_1e8_asy-PE.png"
## Assume 72 pixels/inch and make bitmap twice as large for display resolution
#set terminal png transparent size 2*width*72, 2*height*72 lw 2
#replot

quit
