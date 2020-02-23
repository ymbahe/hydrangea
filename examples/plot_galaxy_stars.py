"""Demonstration script to find particles belonging to a particular galaxy."""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt

# Set script parameters
sim_index = 0               # Which simulation do we want?
first_snapshot_index = 11       # Which snapshot? 29 --> z = 0
second_snapshot_index = 29      # Second snapshot? 29 --> z = 0
ref_snapshot_index = 11
galaxy_id = 10              # Which galaxy to plot?
imsize = 100                 # (Half-)size of analysis region, in kpc
nbins = 100                 # Number of bins for plotting
ptype = 4                   # Look at stars here
plotloc = 'galaxy_remnants.png'  # Where to save the output plot?

# Set up the simulation object for the run we're working with
sim = hy.Simulation(index=sim_index)

# Set up the figure
fig = plt.figure(figsize=(8/0.9, 4))
ax1 = fig.add_axes([0.07, 0.15, 0.4, 0.8])
ax2 = fig.add_axes([0.5, 0.15, 0.4, 0.8])

# Form the (first) files of the required subfind and snapshot catalogues,
# in the first snapshot
subfind_file_ref = sim.get_subfind_file(ref_snapshot_index)

gal_positions = hy.hdf5.read_data(
    sim.gps_loc, 'Centre', read_index=galaxy_id)

shi_ref = hy.hdf5.read_data(sim.fgt_loc, 'SHI',
                            read_index=galaxy_id)[ref_snapshot_index]


def get_log_sigma(snapshot_index):
    """Make a surface density map at a given snapshot."""
    gal_centre = gal_positions[snapshot_index, :]

    # Set up a ReadRegion to extract star particles around the galaxy's
    # position in the current snapshot
    print("\nSetting up <parts> reader:")
    snapshot_file = sim.get_snap_file(snapshot_index)
    parts = hy.ReadRegion(snapshot_file, ptype, gal_centre, imsize,
                          shape='cube')

    # Get particle coordinates relative to subhalo centre (in kpc)
    part_relpos = (parts.Coordinates - gal_centre[None, :]) * 1e3

    # Now filter out only those particles that are actually part of the
    # galaxy in <ref_snapshot_index>:
    ind_in_gal = parts.in_subhalo(shi_ref, subfind_file_ref)

    # Make a simple 2D histogram map of mass surface density
    sigma = (np.histogram2d(part_relpos[ind_in_gal, 1],
                            part_relpos[ind_in_gal, 0], bins=nbins,
                            range=[[-imsize, imsize], [-imsize, imsize]],
                            weights=parts.Mass[ind_in_gal])[0]
             / ((imsize/nbins)**2))

    return np.log10(sigma + 1e-15)


log_sigma_1 = get_log_sigma(first_snapshot_index)
log_sigma_2 = get_log_sigma(second_snapshot_index)

# Plot the surface densities, in log scale
ind_good = np.nonzero(log_sigma_1 >= -14)
vmin, vmax = np.percentile(log_sigma_1[ind_good], [0.01, 99.9])

plt.sca(ax1)
im = plt.imshow(log_sigma_1,
                extent=[-imsize, imsize, -imsize, imsize],
                aspect='equal', interpolation='nearest',
                origin='lower', alpha=1.0, cmap=plt.cm.inferno,
                vmin=vmin, vmax=vmax)

plt.sca(ax2)
im = plt.imshow(log_sigma_2,
                extent=[-imsize, imsize, -imsize, imsize],
                aspect='equal', interpolation='nearest',
                origin='lower', alpha=1.0, cmap=plt.cm.inferno,
                vmin=vmin, vmax=vmax)

# Add a colour bar on the right side
ax2 = fig.add_axes([0.91, 0.15, 0.02, 0.8])
ax2.set_xticks([])
ax2.set_yticks([])
cbar = plt.colorbar(im, cax=ax2, orientation='vertical')
fig.text(0.95, 0.5, r'log$_{10}$ ($\Sigma$ [M$_\odot$ kpc$^{-2}$])',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
