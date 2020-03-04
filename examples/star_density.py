"""Demonstrate loading particles near a given subhalo.

This is an example of how to use the ReadRegion class to quickly read in
particle data within a small sub-volume of the simulation. It creates
a simple stellar surface density map of one galaxy.
"""

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
print("Setting up <subhalo> reader:")
subhalo = hy.SplitFile(subfind_file, 'Subhalo', read_index=sh_index)

# Set up a reading region around subhalo.
# By default, coordinates are expected in Mpc, the same as the default
# output from <subhalo>. From it, we can access all properties of the
# particles in the selected region as attributes of <parts> (analogous to
# <subhalo>). Note that 'size' here refers to the half-side length of
# the cube; with a 'sphere' shape it would be the radius.
print("\nSetting up <parts> reader:")
parts = hy.ReadRegion(snapshot_file,               # (First) file to read
                      ptype,                       # Particle type to read
                      subhalo.CentreOfPotential,   # Centre of region
                      imsize/1e3,                  # Size of region
                      shape='cube')                # Shape of region
print("")

# -------- Aside ----------------------------------------------------------

# With <parts>, we are guaranteed to read in all particles within the
# specified cube, but typically also some particles beyond it will be read.
# We can also specify that we only want particles exactly in the selection
# region, but this takes slightly longer to set up (here using a sphere):
parts_sphere = hy.ReadRegion(snapshot_file, ptype, subhalo.CentreOfPotential,
                             imsize/1e3, shape='sphere', exact=True,
                             verbose=0)  # Silence in the library...

# Now we can, for example, get the total mass or total mass-weighted
# metallicity of all stars within the specified sphere:
galaxy_mstar = parts_sphere.total_in_region('Mass')
galaxy_zstar = parts_sphere.total_in_region('Metallicity', weight_quant='Mass')
print("------------------------------------------------------")
print(f"Subhalo {sh_index} has M_star = {galaxy_mstar:.3e} M_sun, "
      f"Z_star = {galaxy_zstar:.3f}")
print("------------------------------------------------------")

# We can also read one particular data set in a different unit system,
# for example stellar velocities in SI (m/s). Unit conversions are ignored
# for dimensionless quantities:
part_vel_si = parts_sphere.read_data("Velocity", units="SI")
part_ids = parts_sphere.read_data("ParticleIDs", units="cgs")
print("----------------------------------------------------------------------")
print(f"Particle ID={part_ids[0]} has velocity "
      f"{part_vel_si[0, 0]:.3f}/{part_vel_si[0, 1]:.3f}/"
      f"{part_vel_si[0, 2]:.3f} m/s")
print("----------------------------------------------------------------------")

# ------- Back to main script purpose -------------------------------------

# Get particle coordinates relative to subhalo centre (in kpc)
print("\nNote how <parts> now reads the 'Coordinates' and 'Mass' data sets:")
part_relpos = (parts.Coordinates - subhalo.CentreOfPotential[None, :]) * 1e3

# Make a simple 2D histogram map of mass surface density
sigma = np.histogram2d(part_relpos[:, 1], part_relpos[:, 0], bins=nbins,
                       range=[[-imsize, imsize], [-imsize, imsize]],
                       weights=parts.Mass)[0] / ((imsize/nbins)**2)

# Plot the histogram, in log scale
fig = plt.figure(figsize=(4/0.8, 4))
ax1 = fig.add_axes([0.15, 0.15, 0.65, 0.8])
log_sigma = np.log10(sigma + 1e-5)  # Add small bias to avoid NaNs
ind_good = np.nonzero(log_sigma >= -4)
vmin, vmax = np.percentile(log_sigma[ind_good], [0.01, 99.9])
im = plt.imshow(log_sigma,
                extent=[-imsize, imsize, -imsize, imsize],
                aspect='equal', interpolation='nearest',
                origin='lower', alpha=1.0, cmap=plt.cm.inferno,
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
