reset
#
set terminal pslatex color rotate auxfile size 3.60,4.80 8
#
Folder = '../Instability/Output/'
#
DataFile = 'Output_Biaxial_Comp_k16_Dcircls_49_a-57_Plot_Gen.txt'
#
unset key
#set grid back
#
# Line styles of plotting DEM results
set style line 1 dt 1 lw 2.0 lc rgb "#8B0000" # dark red       Spheres
set style line 2 dt 12 lw 2.5 lc rgb "#191970" # midnight blue  Ovoids_0800
set style line 3 dt 3 lw 2.0 lc rgb "#006400" # dark green     Ovoids_0625
set style line 4 dt 4 lw 2.0 lc rgb "black"   #                Ovoids_0500
#
set style line 5 dt 5 lw 2.5 lc rgb "purple"
set style line 6 dt 1 lw 3 lc rgb "blue"
set style line 7 dt 5 lw 2.0 lc rgb "orangered4"
#
# black lines
set style line 8 dt 2  lw 1 lc rgb "black"    # Thin black lines
set style line 9 dt 3  lw 1.5 lc rgb "black"    # Zero axis  lines
set style line 10 dt 2 lw 1.4 lc rgb "black"    # Zero axis  lines
set style line 11 dt 6 lw 2.0 lc rgb "#191970"
set style line 12 dt 1 lw 1.4 lc rgb "green"
#
# Experimental (Tsukamoto) lines
set style line 20 dt 3 lw 1.0 lc rgb "black"
set style line 21 dt 3 lw 1.4 lc rgb "green"
set style line 22 dt 1 lw 1.0 lc rgb "black"
set style line 23 dt 1 lw 1.0 lc rgb "black"
set style line 24 dt 2 lw 0.5 lc rgb "black"
#
set border lw 2.0
#
# set logscale x
# set xzeroaxis ls 20
# set xtics format ""
#
set multiplot
set lmargin 10.0
set bmargin 0.5
set grid xtics ytics mxtics ls 24
#
# Set labels here
# set label '\parbox{3.00in}{\centering Spectral acceleration, g}' \
#     at screen 0.05,0.5 center front rotate
#
set size 0.8, 0.25
set origin 0.0, 0.75
#
set pointsize 0.2
#
unset xlabel
#set xlabel 'Compressive strain, $-\sigma_{11}$' offset 0,0.3
set ylabel '\parbox{1.00in}{\centering Stress ratio,\\ $q/p$}' offset -1.3,0
#
set xrange[0:4.4];  set xtics offset 0,0.0 2; set mxtics 2; set format x ""
set yrange[0:1.25]; set ytics offset 0.3,0 0.5; set mytics 5; set format y "%.1f"
#
# set label '$2Nb/a$' at 2.2,0.86
#set key at 41.7,-0.8 left top Left reverse samplen 3.5 spacing 1.20 nobox
#
set arrow from 1,screen 0.480 to 1,screen -0.04 nohead ls 24
set arrow from 2,screen 0.480 to 2,screen -0.04 nohead ls 24
set arrow from 3,screen 0.480 to 3,screen -0.04 nohead ls 24
set arrow from 4,screen 0.480 to 4,screen -0.04 nohead ls 24
#
plot Folder.DataFile using (100*$3):4 title 'q/p ratio' with points ls 1 pt 7
#
unset arrow
#
set size 0.8,0.19
set origin 0, 0.55
#
#set format x "%.2f"
set format x "%g"
#
set yrange[-0.01:0.01]; set ytics offset 0.3,0 0.01; set mytics 2; set format y "%g"
#
set xlabel 'Compressive strain, $-\varepsilon_{11}$, percent' offset 0,0.4
set ylabel '\parbox{1.00in}{\centering Volume strain,\\ $dV/V_{\text{o}}$}' offset 0.5,0
#
set xzeroaxis ls 23
#
plot Folder.DataFile using (100*$3):5 title 'dv' with points ls 1 pt 7
#
set noxzeroaxis
unset xlabel
unset ylabel; unset xtics; unset mxtics
unset ylabel; unset ytics; unset mytics
#
set border lw 1.0
#
set yrange[0:1]
#
# non-P-matrix
set size 0.8,0.07; set origin 0, 0.42
set pointsize 0.2
#
plot Folder.DataFile using (100*$3):($11==0 ? 0.9*rand(0)+0.05 : 1/0) title 'non-P-matrix' with points ls 1 pt 7
#
# non-R-matrix
set size 0.8,0.07; set origin 0, 0.37
set pointsize 0.2
#
plot Folder.DataFile using (100*$3):($22==0 ? 0.9*rand(0)+0.05 : 1/0) title 'non-R-matrix' with points ls 1 pt 7
#
# Degenerate matrix
set size 0.8,0.07; set origin 0, 0.32
set pointsize 0.2
#
plot Folder.DataFile using (100*$3):($13==1 ? 0.9*rand(0)+0.05 : 1/0) title 'degenerate' with points ls 1 pt 7
#
# Type I bifurcation
set size 0.8,0.07; set origin 0, 0.26
set pointsize 0.3
#
plot Folder.DataFile using (100*$3):($14>1 || ($14>0 && $15>0) ? 0.25*rand(0)+0.375 : 1/0) \
     title 'Type I' with points ls 1 pt 7
#
# Type II bifurcation
set size 0.8,0.07; set origin 0, 0.21
set pointsize 0.3
#
plot Folder.DataFile using (100*$3):($15>0 ? 0.30*rand(0)+0.35 : 1/0) \
     title 'Type II' with points ls 1 pt 7
#
# Path instability
set size 0.8,0.07; set origin 0, 0.16
set pointsize 0.35
#
plot Folder.DataFile using (100*$3):(($14)>1 ? 0.30*rand(0)+0.35 : 1/0) \
     title 'Path instability' with points ls 1 pt 7, \
     Folder.DataFile using (100*$3):(($15)>=1 ? 0.30*rand(0)+0.35 : 1/0) \
     title 'Path instability' with points ls 6 pt 7
#
# Path-sensitive
set size 0.8,0.07; set origin 0, 0.11
set pointsize 0.3
#
plot Folder.DataFile using (100*$3):($17>0 ? 0.30*rand(0)+0.35 : 1/0) \
     title 'Path-sensitive' with points ls 1 pt 7
#
# Path-incongruity
set size 0.8,0.07; set origin 0, 0.06
set pointsize 0.2
#
plot Folder.DataFile using (100*$3):($19>0 ? 0.50*rand(0)+0.25 : 1/0) \
     title 'Path-incongruity' with points ls 1 pt 7
#
# Neutral stability
set size 0.8,0.07; set origin 0, 0.01
set pointsize 0.3
#
plot Folder.DataFile using (100*$3):($18>0 ? 0.30*rand(0)+0.35 : 1/0) \
     title 'Neutral stability' with points ls 1 pt 7
#
# Instability of equilibrium
set size 0.8,0.07; set origin 0, -0.04
set pointsize 0.3
#
plot Folder.DataFile using (100*$3):($20>0 ? 0.9*rand(0)+0.05 : 1/0) \
     title 'Instability of equilibrium' with points ls 1 pt 7
#
set arrow from screen 0,0 to screen 1,1 nohead ls 1
#
unset multiplot
