import numpy as np
import matplotlib.pyplot as plt
import horpy as hp

# -----------------------------
# PURPOSE
# -----------------------------

# With custom grid spacings, static mesh coarsening, etc,
# it is difficult to picture the grid structure in your head.
# This script generates a direct visualization of the grid
# structure for spherical coordinate systems.
# Only for x-z plane.

# -----------------------------
# USER INPUTS
# -----------------------------

# Specify size of plotting region (Rsun)
r_max = 2.

# Set up HORMONE variables as it would in the simulation
# Give the path to the work directory as the argument
hp.readbin_python.setup_python('../')

# You can read grid information by feeding in the gridfile
hp.readbin_python.read_gridfile('gridfile.bin')

ie = hp.grid.ie
je = hp.grid.je
rsun=hp.constants.rsun
xi1 = hp.readbin_python.get_interior_1d(hp.grid.xi1)/rsun
xi2 = hp.grid.xi2[1:je+3]

# FMR layers: number of radial cells per layer
fmr_max = hp.grid.fmr_max
fmr_lvl = hp.grid.fmr_lvl[1:fmr_max+1]

# -----------------------------
# BUILD LAYER STRUCTURE
# -----------------------------

# Angular resolution per layer
polar_cells_per_layer = [
    je // (2 ** (fmr_max - k - 0)) for k in range(fmr_max)
]

# Radial layers
radial_layers = []
i0 = 0
for ncell in fmr_lvl:
    i1 = i0 + ncell
    radial_layers.append((i0, i1))
    i0 = i1

# Add outer full-resolution region
if i0 < ie:
    radial_layers.append((i0, ie))
    polar_cells_per_layer.append(je)

# -----------------------------
# PLOTTING
# -----------------------------

fig, ax = plt.subplots(figsize=(6, 6))

# Dense theta grid for smooth arcs
theta_dense = np.linspace(0.0, 0.5 * np.pi, 400)

for (i_start, i_end), npol in zip(radial_layers, polar_cells_per_layer):

    # Polar interfaces for this layer
    stride = je // npol
    theta_interfaces = xi2[::stride]
    if theta_interfaces[-1] != xi2[-1]:
        theta_interfaces = np.append(theta_interfaces, xi2[-1])

    # Radial interfaces
    r_interfaces = xi1[i_start : i_end + 1]

    # ---- constant-theta lines (radial rays) ----
    for theta in theta_interfaces:
        r = r_interfaces
        R = r * np.sin(theta)
        Z = r * np.cos(theta)
        ax.plot(R, Z, color="black", lw=0.4)

    # ---- constant-r lines (circular arcs) ----
    for r in r_interfaces:
        R = r * np.sin(theta_dense)
        Z = r * np.cos(theta_dense)
        ax.plot(R, Z, color="black", lw=0.4)

# -----------------------------
# FINAL TOUCHES
# -----------------------------

ax.set_aspect("equal")
ax.set_xlim(0, r_max)
ax.set_ylim(0, r_max)

ax.set_xlabel(r"Radius ($R_\odot$)")
ax.set_ylabel("Polar axis")

ax.tick_params(direction="in")
plt.tight_layout()
plt.show()

