"""Demonstration script to trace particles between two snapshots."""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt
import sys
from pdb import set_trace

# Set script parameters
sim_index = 0                   # Which simulation do we want?
snap_index_ref = 29             # Snapshot for particle selection (z = 0)
snap_index_plot = 11            # Which snapshot to analyse? (z ~ 1)
sel_size = 0.2                  # Selection radius of particles, in r500c
temp_range = [2.5, 9.0]         # Plot range in (log) temperature [K]
nH_range = [-9.0, 3.0]          # Plot range in (log) nH [cm^-3]
scale_range = [9.0, 14.5]       # Scale range of plot
nbins = 100                     # Number of bins per axis
plotloc = 'cluster_origin.png'  # Where to save the output plot?

# Set up the simulation object for the run we're working with
sim = hy.Simulation(index=sim_index)

# Form the (first) files of the required subfind and snapshot catalogues
subfind_file_ref = sim.get_subfind_file(snap_index_ref)
snapshot_file_ref = sim.get_snap_file(snap_index_ref)
snapshot_file_plot = sim.get_snap_file(snap_index_plot)

# Set up a reader for the subhalo catalogue, specifying that we are only
# interested in one particular object (to speed things up). From this,
# we can access all properties of the selected subhalo as attributes
# of <subhalo>; they are loaded when first required. By default, quantities
# are returned in 'astronomically sensible' units (e.g. Mpc, Gyr, km/s).
fof_cl = hy.SplitFile(subfind_file_ref, 'FOF', read_index=0)

# Set up a reader to select gas particles within the desired aperture
# around the cluster centre at z = 0
gas_ref = hy.ReadRegion(snapshot_file_ref, 0, fof_cl.GroupCentreOfPotential,
                        sel_size * fof_cl.Group_R_Crit500, exact=True)

# Set up a second reader, to read the gas particle catalogue in the
# second (plotting) snapshot. We are interested in particles that were
# in the cluster centre at z = 0, but we don't know where they were in
# this snapshot, so we have to read in the entire catalogue
# (in principle, we could also try successively larger apertures around
# the cluster progenitor, see plot_thermal_history.py example)
gas_plot = hy.SplitFile(snapshot_file_plot, part_type=0)

# Now we need to identify the particles in the plot snapshot that were
# in the selected region in the ref snapshot. Here, we use the
# find_id_indices functtion; we could also have used the Gate class
# as in the snipshot_stellar_age.py example.

print(f"\nMatching particles between snaps {snap_index_ref} and "
      f"{snap_index_plot}, also reading ParticleIDs...")
ind_in_plot, ind_matched = hy.crossref.find_id_indices(gas_ref.ParticleIDs,
                                                       gas_plot.ParticleIDs)

# ind_matched lists all particles that could be matched. This should be
# all, but let's check to be sure:
if len(ind_matched) != len(ind_in_plot):
    print(f"Problem: could only match {len(ind_matched)} out of "
          f"{len(ind_in_plot)} particles. Please investigate.")
    set_trace()

# Now plot a phase diagram of the particles, in analogy to the
# plot_phase_diagram.py example

# Prepare the plot
fig = plt.figure(figsize=(5/0.8, 4))
ax1 = fig.add_axes([0.12, 0.15, 0.65, 0.8])

# Create the 2D histogram. The particle temperatures, densities, and
# masses are read in implicitly as they are being accessed.
histogram, xe, ye = np.histogram2d(
    np.log10(gas_plot.Temperature[ind_in_plot]),
    np.log10(gas_plot.Density[ind_in_plot]),
    weights=gas_plot.Mass[ind_in_plot],
    bins=[nbins, nbins],
    range=[temp_range, nH_range])

# Avoid NaN values in log by adding a *very* small offset to all values
histogram += sys.float_info.min

# Plot the 2D histogram
da = (xe[1]-xe[0]) * (ye[1]-ye[0])
im = plt.imshow(np.log10(histogram/da), cmap=plt.cm.inferno,
                extent=[*nH_range, *temp_range], origin='lower', aspect='auto',
                vmin=scale_range[0], vmax=scale_range[1])

ax1.set_xlim(nH_range)
ax1.set_ylim(temp_range)
ax1.set_xlabel(r'$\log_{10}\,(n_\mathrm{H}\,[\mathrm{cm}^{-3}])$')
ax1.set_ylabel(r'$\log_{10}\,(T\,[\mathrm{K}])$')

# Add a colour bar on the right side
ax2 = fig.add_axes([0.81, 0.15, 0.05, 0.8])
ax2.set_xticks([])
ax2.set_yticks([])
cbar = plt.colorbar(im, cax=ax2, orientation='vertical')
fig.text(0.94, 0.55, r'log$_{10}$ (d$M_\mathrm{gas}$ / '
         r'(dlog$_{10}\,T$ dlog$_{10}\,n_\mathrm{H}$) [M$_\odot$])',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
