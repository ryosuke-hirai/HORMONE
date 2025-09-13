import horpy as hp

hp.readbin_python.setup_python('../')
hp.readbin_python.read_gridfile('gridfile.bin')
hp.readbin_python.read_binfile('bin00000000000ms.dat')

darray =hp.readbin_python.d.copy()
parray =hp.readbin_python.p.copy()
x1array=hp.readbin_python.x1.copy()
x2array=hp.readbin_python.x2.copy()
x3array=hp.readbin_python.x3.copy()
print(darray)
#exec(open('plot.py').read())
