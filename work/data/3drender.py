#!../../pyenv/bin/python
import horpy as hp
import numpy as np
import pyvista as pv
import glob
import os
import re

# ==========================
# Step 1: Setup
# ==========================
hp.readbin_python.setup_python('../')
hp.readbin_python.read_gridfile('gridfile.bin')

# Collect all bin files
binfiles = sorted(glob.glob("bin*min.dat"))
print(f"Found {len(binfiles)} bin files total")

binfiles = binfiles[0:10:1]
print(f"Processing {len(binfiles)} selected files")

# Output directory
os.makedirs("frames", exist_ok=True)

# Fixed camera position for consistency
camera_pos = [(-7e12, -7e12, 7e12), (0, 0, 0), (0, 0, 1)]

# Some physical constants
clight = 2.99792458e10
kbol = 1.380649e-16
hplanck = 6.62607015e-27
sigma=2*np.pi**5*kbol**4/(15*clight**2*hplanck**3)
arad = 4*sigma/clight
Navo = 6.02214076e23
Rgas = kbol*Navo

# ==========================
# Step 2: Loop over bin files
# ==========================
for binfile in binfiles:
#    print(f"Processing {binfile} ({i+1}/{len(binfiles)})")
    print(f"Processing {binfile} ...")

    # Extract the number part (e.g., "00000144000" from "bin00000144000min.dat")
    match = re.search(r'bin(\d+)min\.dat', binfile)
    if match:
        filenum = match.group(1)
    else:
        filenum = "unknown"

    hp.readbin_python.read_binfile(binfile)

    # Extract interior arrays
    darray  = hp.readbin_python.get_interior_3d(hp.physval.d)
    Tarray  = hp.readbin_python.get_interior_3d(hp.physval.t)
    parray  = hp.readbin_python.get_interior_3d(hp.physval.p)
    x1array = hp.readbin_python.get_interior_1d(hp.grid.x1)   # r
    x2array = hp.readbin_python.get_interior_1d(hp.grid.x2)   # theta
    x3array = hp.readbin_python.get_interior_1d(hp.grid.x3)   # phi
    sink_x  = hp.readbin_python.sink_x.copy()

    # ==========================
    # Step 3: Downsample safely
    # ==========================
    r_slice   = slice(None, None, 8)
    th_slice  = slice(None, None, 4)
    phi_slice = slice(None, None, 4)

    r   = x1array[r_slice]
    th  = x2array[th_slice]
    phi = x3array[phi_slice]
    rho = darray[r_slice, th_slice, phi_slice]
    temp= Tarray[r_slice, th_slice, phi_slice]
    pres= parray[r_slice, th_slice, phi_slice]
    ent = 1/0.62*np.log(temp**1.5/rho) + 4*arad*temp**3/(3*rho)

    plotval = ent

    # ==========================
    # Step 4: Extend grid
    # ==========================
    phi_ext = np.append(phi, phi[0])
    val_phiext = np.concatenate([plotval, plotval[:, :, :1]], axis=2)

    th_ext = np.append(0.0, th)
    val0 = np.expand_dims(val_phiext[:, 0, :], axis=1)
    val_thext = np.concatenate([val0, val_phiext], axis=1)

    th_mirror = np.pi - th_ext
    th_full = np.concatenate([th_ext, th_mirror[::-1]])
    val_mirror = val_thext[:, ::-1, :]
    val_full = np.concatenate([val_thext, val_mirror], axis=1)

    val_log = np.log10(np.clip(val_full, 1e-30, None))

    # ==========================
    # Step 5: Build StructuredGrid and visualize
    # ==========================
    rr, tt, pp = np.meshgrid(r, th_full, phi_ext, indexing='ij')
    x = rr * np.sin(tt) * np.cos(pp)
    y = rr * np.sin(tt) * np.sin(pp)
    z = rr * np.cos(tt)

    grid = pv.StructuredGrid(x, y, z)
    grid['log_density'] = val_log.ravel(order='F')

#    vmin, vmax = -9, -6
    vmin, vmax = np.percentile(val_log, [5, 95])
    iso_values = np.linspace(vmin, vmax, 5)
#    iso_values = np.append(iso_values, -4)
    contours = grid.contour(isosurfaces=iso_values, scalars='log_density')

    sink_coord = sink_x[:, 1]  # second sink
    point_mesh = pv.PolyData(sink_coord.reshape(1, 3))

    p = pv.Plotter(off_screen=True)

    opacity = [1, 0.3, 0.05, 0.02, 0.01]
    p.add_mesh(contours, cmap='plasma', opacity=0.1)
    p.add_mesh(point_mesh, color='black', point_size=10, render_points_as_spheres=True, opacity=1)
    
    p.camera_position = camera_pos
#    p.set_background("black")

    # Output filename (denXXXXX.png)
    outname = f"frames/den{filenum}min.png"

    # Render and save
#    p.show()
    p.show(auto_close=False)
    p.screenshot(outname, window_size=[960, 560])
    p.close()

    print(f"Saved {outname}")
