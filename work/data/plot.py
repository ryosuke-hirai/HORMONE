# Import horpy.so
import horpy as hp

# Set up HORMONE variables as it would in the simulation
# Give the path to the work directory as the argument
hp.readbin_python.setup_python('../')

# You can read grid information by feeding in the gridfile
hp.readbin_python.read_gridfile('gridfile.bin')

# You can read binary dumps by feeding in the binXX.dat files
hp.readbin_python.read_binfile('bin00000000000ms.dat')

# You can directly access all the HORMONE variables in src/main/modules.f90
# Most of the relevant variables should be in the following modules
# grid: Grid-related quantities
# physval: Physical quantities
# gravmod: Gravity-related quantities
darray_full = hp.physval.d.copy()
x1array_full = hp.grid.x1.copy()
#grvphi_full = hp.gravmod.grvphi.copy()

# There are helper functions to chop off ghost cells
# get_interior_1d: For 1D variables like x1, x2, x3, etc
# get_interior_3d: For 3D variables like d, p, e, v1, v2, v3, etc
# get_interior_4d: For 4D variables like spc
darray  = hp.readbin_python.get_interior_3d(hp.physval.d)
parray  = hp.readbin_python.get_interior_3d(hp.physval.p)
x1array = hp.readbin_python.get_interior_1d(hp.grid.x1)
x2array = hp.readbin_python.get_interior_1d(hp.grid.x2)
x3array = hp.readbin_python.get_interior_1d(hp.grid.x3)
#spcarray = hp.readbin_python.get_interior_4d(hp.physval.spc)

# Sink particle information is stored directly under readbin_python
# sink_x(3,n): 2D array where first index is for x1,x2,x3 coordinates and second index is particle number
# sink_v(3,n): 2D array where first index is for v1,v2,v3 velocity elements and second index is particle number
# sink_mass(n): 1D array for sink particle masses
# sink_mdot(n): 1D array for sink particle mass accretion rates
#sink_x_array = hp.readbin_python.sink_x.copy()
#sink_v_array = hp.readbin_python.sink_v.copy()
#sink_mass_array = hp.readbin_python.sink_mass
#sink_mdot_array = hp.readbin_python.sink_mdot

#exec(open('plot.py').read())
