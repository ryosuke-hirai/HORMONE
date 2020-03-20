reset
clear

set view map
set size ratio 0.5
set yr [0:*]

amp = 200
spl 'data/plt09505.dat' ev ::2::129 u (-$4*cos($5)):($4*sin($5)):7 w pm3d , '' ev 3 u (-$4*cos($5)):($4*sin($5)):(0):((-cos($5)*$10+sin($5)*$11)*amp):((cos($5)*$11+sin($5)*$10)*amp):(0) w vector lc 5