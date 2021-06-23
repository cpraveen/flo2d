       reset
       set term svg size 800,800
       set out 'flo.svg'

       set multiplot

       set size 0.5,0.5
       set origin 0.0,0.5
       set xran[-11:11]
       set noxtics
       set nokey
#      set bmargin 1
       set title "Pressure/Friction coeffient"
#      set y2tics
       set ylabel '-Cp'
#      set y2label 'Cf'
       plot 'WALL.DAT' u 1:3 w p lw 1 pt 6

       reset
       set size 0.5,0.5
       set origin 0.5,0.5
       set auto
       set noxtics
       set nokey
       set logscale y
       set title "Residue"
       plot 'FLO.RES' u 1:2 w l
       set nologscale y

       reset
       set size 0.75,0.75
       set origin 0.2,0.0
       set xran[-11:11]
       set yran[-3:3]
       set size ratio -1
       set noxtics
       set noytics
       set nokey
#      set title "Mach number"
#      set bmargin 1
       plot 'FLO.M' w l,'BD.DAT' w l lt 1 lw 2

       reset
       set size 0.75,0.75
       set origin 0.2,-0.2
       set xran[-11:11]
       set yran[-3.0:3.0]
       set size ratio -1
       set noxtics
       set noytics
       set nokey
#      set bmargin 1
#      set title "Pressure"
       plot 'FLO.P' w l,'BD.DAT' w l lt 1 lw 2

       unset multiplot
#      set term x11
