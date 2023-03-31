reset
#clear

box = 1
set xr [0:box]
set yr [-box:box]

pltfile(i)=sprintf('<paste gridfile.dat plt%11.11dms.dat',i)

xaxis='x1'

do for [i=200:200:10]{

timelabel=sprintf('%d ms',i)
set label 1 timelabel at graph 0.5,graph 0.95 center front

    plot pltfile(i) u ($2):(column('d'))  w l t 'd',\
	 pltfile(i) u ($2):(column('p'))  w l t 'p',\
	 pltfile(i) u ($2):(column('b1')) w l t 'b1',\
	 pltfile(i) u ($2):(column('b2')) w l t 'b2',\
	 pltfile(i) u ($2):(column('b3')) w l t 'b3',\
	 pltfile(i) u ($2):(column('v1')) w l t 'v1',\
	 pltfile(i) u ($2):(column('v2')) w l t 'v2',\
	 pltfile(i) u ($2):(column('v3')) w l t 'v3'

	 

}
