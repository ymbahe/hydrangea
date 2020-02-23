"""Demonstration script to find particles belonging to a particular galaxy."""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt
from pdb import set_trace

# Set script parameters
sim_index = 0                    # Which simulation do we want?
first_snapshot_index = 8         # Which snapshot? 29 --> z = 0
second_snapshot_index = 12       # Second snapshot? 29 --> z = 0
ref_snapshot_index = 8
galaxy_id = 1808                   # Which galaxy to plot?
imsize = 100                     # (Half-)size of analysis region, in kpc
nbins = 100                      # Number of bins for plotting
ptype = 4                        # Look at stars here
plot_range = [6.0, 9.5]
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


def plot_log_sigma(snapshot_index, ax):
    """Plot a surface density map at a given snapshot."""
    gal_centre = gal_positions[snapshot_index, :]

    # Set up a ReadRegion to extract star particles around the galaxy's
    # position in the current snapshot
    print("\nSetting up <parts> reader:")

    snapshot_file = sim.get_snap_file(snapshot_index)
    parts = hy.ReadRegion(snapshot_file, ptype, gal_centre, imsize/1e3,
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

    plt.sca(ax)
    im = plt.imshow(np.log10(sigma + 1e-15),
                    extent=[-imsize, imsize, -imsize, imsize],
                    aspect='equal', interpolation='nearest',
                    origin='lower', alpha=1.0, cmap=plt.cm.inferno,
                    vmin=plot_range[0], vmax=plot_range[1])

    return im


im = plot_log_sigma(first_snapshot_index, ax1)
im = plot_log_sigma(second_snapshot_index, ax2)


# Add a colour bar on the right side
ax3 = fig.add_axes([0.91, 0.15, 0.02, 0.8])
ax3.set_xticks([])
ax3.set_yticks([])
cbar = plt.colorbar(im, cax=ax3, orientation='vertical')
fig.text(0.95, 0.5, r'log$_{10}$ ($\Sigma$ [M$_\odot$ kpc$^{-2}$])',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
