if (!exists("PREFIX")) PREFIX = "phase_separation_physical"
if (!exists("OUT"))    OUT    = "phase_separation_physical_psi.gif"
if (!exists("FIRST"))  FIRST  = 1000
if (!exists("LAST"))   LAST   = 200000
if (!exists("EVERY"))  EVERY  = 1000
if (!exists("DT"))     DT     = 5.0e-3

if (!exists("NX")) NX = 128
if (!exists("NY")) NY = 128

# For a=-3, u=1, the homogeneous minima are psi = +/- sqrt(3).
if (!exists("CBMIN")) CBMIN = -1
if (!exists("CBMAX")) CBMAX =  1

set datafile commentschars "#"
set terminal gif animate optimize delay 8 size 900,820
set output OUT

set size ratio -1
set xrange [-0.5:NX-0.5]
set yrange [-0.5:NY-0.5]
# set cbrange [CBMIN:CBMAX]

set xlabel "x grid"
set ylabel "y grid"
set cblabel "psi"
set title ""

set palette defined ( \
    CBMIN "#313695", \
    0.0   "#ffffbf", \
    CBMAX "#a50026" )

unset key
set view map

do for [step=FIRST:LAST:EVERY] {
    file = sprintf("%s_step%07d.snapshot", PREFIX, step)
    t = step * DT
    set title sprintf("phi^4 phase separation: psi(x,y), step = %d, t = %.3g", step, t)
    plot file using 1:2:7 with image
}

set output
