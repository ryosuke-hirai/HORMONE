reset
#clear

box = 1
set xr [0:box]
#set yr [-box:box]

pltfile(i)=sprintf('<paste gridfile.dat plt%11.11dms.dat',i)

xaxis='x1'

do for [i=300:300:10]{

timelabel=sprintf('%d ms',i)
set label 1 timelabel at graph 0.5,graph 0.95 center front

    plot pltfile(i) u ($2):(column('d'))  w l t 'd',\
	 pltfile(i) u ($2):(column('p'))  w l t 'p',\
	 pltfile(i) u ($2):(column('v1')) w l t 'v1'

	 

}
