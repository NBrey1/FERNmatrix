set macro  # Enable macro definition

set termoption dashed   # May be required for dashed .png

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
width = 8.5
height = 8.5

# x-axis resolution
set samples 1000

# Line styles.  
# For lines: plot x with lines ls 1
# For points: plot x with points ls 1

set style line 1 lc rgb 'black' pt 5 lw 1.0   # fill square
set style line 2 lc rgb myred pt 6 lw 1.0  # open circle
set style line 3 lc rgb 'blue' pt 9 lw 1.0  # fill triangle
set style line 4 lc rgb mygreen pt 6  lw 1.0 # open circle
set style line 5 lc rgb mypurple pt 6 lw 1.25 # open circle
set style line 6 lc rgb myviolet  pt 7   # fill circle
set style line 7 lc rgb mybrown  pt 7   # fill circle
set style line 8 lc rgb myorange  pt 7   # fill circle
set style line 9 lc rgb myblue pt 7   # fill circle
set style line 10 lc rgb 'green' pt 6   # open circle
set style line 11 lc rgb 'red' pt 7   # open square
set style line 12 lc rgb 'gray' pt 4   # open square
set style line 13 lc rgb 'gold' pt 4   # open square
set style line 14 lc rgb 'orange' pt 4   # open square
set style line 15 lc rgb 'magenta' pt 4   # open square

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'log10 Time (s)' textcolor rgb tic_color #font "Arial,32"
set ylabel 'log10 Value' textcolor rgb tic_color #font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

# Use scientific exponential notation on axes

#set format x "%3.1t{x}10^{%L}"  # scientific
#set format y "%3.1t{x}10^{%L}"  # scientific

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="Asy T9=7 rho=1e8 (no PF; RG 9)"
ds = ds.""
set title ds textcolor rgb title_color #font "Arial,22"


# -------- Axis ranges and ticmarks -----------

xlow = -8.5
xup = -7.6
xtics = 0.1   # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -3.0#-2.87 #5.8
yup = -2.55#-2.7 #6.2
ytics = 0.05      # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic
set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------

set title ds textcolor rgb title_color

file1 = "gdb-C++_Asy_36Ar.out"

# Edit the following plot commands to correspond to data
# read in from data file and convert quantities to
# log10 of quantities to make log-log plot

#plot file1 using (log10($6)):(log10($7)) with lines ls 1 lw 1.0 dashtype 1 title "Asy Y0[36Ar]"
plot file1 using (log10($6)):(1.0*log10($8)) with lines ls 2 lw 1.0 dashtype 1 title "log10(Asy Y[36Ar])"
replot file1 using (log10($6)):(log10($9)-8.8) with lines ls 3 lw 1.0 dashtype 1 title "log10(Asy flux[31])-8.8"

file2 = "gdb-C++_PE_36Ar.out"

#replot file2 using (log10($6)):(log10($7)) with lines  ls 1 lw 2.0 dashtype 2 title "PE Y0[36Ar]"
replot file2 using (log10($6)):(1.0*log10($8)) with lines  ls 2 lw 1.0 dashtype 2 title "log10(PE Y[36Ar])"
replot file2 using (log10($6)):(log10($9)-8.8) with lines ls 3 lw 1.0 dashtype 2 title "log10(PE flux[31])-8.8"


# Reset font sizes for .eps and .png output2

set key top right font "Arial,22"
set title ds textcolor rgb title_color font "Arial,22"
set xlabel 'log10 Time (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'log10 Value' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "diagnostics_36Ar_log.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file

# Plot to PNG file

#set out "diagnostics_36Ar_log.png"
## Assume 72 pixels/inch and make bitmap twice as large for display resolution
#set terminal png transparent size 2*width*72, 2*height*72 lw 2
#replot

quit
