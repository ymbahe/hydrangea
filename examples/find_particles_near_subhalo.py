"""Demonstrate loading particles near a given subhalo."""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt

# Set script parameters
sim_index = 0               # Which simulation do we want?
snap_index = 29             # Which snapshot? 29 --> z = 0
sh_index = 8                # Which subhalo? 
imsize = 50                 # (Half-)size of analysis region, in kpc
nbins = 100                 # Number of bins for plotting
ptype = 4                   # Look at stars here
plotloc = 'galaxy_map.png'  # Where to save the output plot?

# Set up the simulation object for the run we're working with
sim = hy.Simulation(index=sim_index)

# Form the (first) files of the required subfind and snapshot catalogues
subfind_file = sim.get_subfind_file(snap_index)
snapshot_file = sim.get_snap_file(snap_index)

# Set up a reader for the subhalo catalogue, specifying that we are only
# interested in one particular object (to speed things up). From this,
# we can access all properties of the selected subhalo as attributes
# of <subhalo>; they are loaded when first required. By default, quantities
# are returned in 'astronomically sensible' units (e.g. Mpc, Gyr, km/s). 
subhalo = hy.SplitFile(subfind_file, 'Subhalo', read_index=sh_index)

# Set up a reading region around subhalo.
# By default, coordinates are expected in Mpc, the same as the default
# output from <subhalo>. From it, we can access all properties of the
# particles in the selected region as attributes of <parts> (analogous to
# <subhalo>).
parts = hy.ReadRegion(snapshot_file,               # (First) file to read
                      ptype,                       # Particle type to read
                      subhalo.CentreOfPotential,   # Centre of region
                      imsize/1e3,                  # Size of region
                      shape='cube')                # Shape of region


# Get particle coordinates relative to subhalo centre (in kpc)
part_relpos = (parts.Coordinates - subhalo.CentreOfPotential[None, :]) * 1e3

# Make a simple 2D histogram map of mass surface density
sigma = np.histogram2d(part_relpos[:, 1], part_relpos[:, 0], bins=nbins,
                       range=[[-imsize, imsize], [-imsize, imsize]],
                       weights=parts.Mass)[0] / (imsize/nbins**2)

# Plot the histogram, in log scale
fig = plt.figure(figsize=(4/0.8, 4))
ax1 = fig.add_axes([0.15, 0.15, 0.65, 0.8])
log_sigma = np.log10(sigma + 1e-5)  # Add small bias to avoid NaNs
ind_good = np.nonzero(log_sigma >= -4)
vmin, vmax = np.percentile(log_sigma[ind_good], [0.01, 99.9])
im = plt.imshow(log_sigma,
                extent = [-imsize, imsize, -imsize, imsize],
                aspect = 'equal', interpolation = 'nearest', 
                origin = 'lower', alpha = 1.0, cmap=plt.cm.inferno,
                vmin=vmin, vmax=vmax)

ax1.set_xlim((-imsize, imsize))
ax1.set_ylim((-imsize, imsize))
ax1.set_xlabel(r'$\Delta x$ [kpc]')
ax1.set_ylabel(r'$\Delta y$ [kpc]')

# Add a colour bar on the right side
ax2 = fig.add_axes([0.81, 0.15, 0.05, 0.8])
ax2.set_xticks([])
ax2.set_yticks([])
cbar = plt.colorbar(im, cax=ax2, orientation='vertical')
fig.text(0.95, 0.5, r'log$_{10}$ ($\Sigma$ [M$_\odot$ kpc$^{-2}$])',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False)




