"""Demonstration script to plot the phase diagram of one simulation.

It illustrates using the SplitFile class to read in an entire particle
catalogue for a snapshot.
"""

import hydrangea as hy
import numpy as np
import matplotlib.pyplot as plt
import sys

# Set script parameters
sim_index = 14                 # Which simulation do we want?
snap_index = 29                # Which snapshot do we want to analyse?
temp_range = [2.5, 9.0]        # Plot range in (log) temperature [K]
nH_range = [-9.0, 3.0]         # Plot range in (log) nH [cm^-3]
scale_range = [9.0, 14.5]      # Scale range of plot
nbins = 100                    # Number of bins per axis
plotloc = 'phase_diagram.png'  # Where to save the output plot?

# Prepare the plot
fig = plt.figure(figsize=(5/0.8, 4))
ax1 = fig.add_axes([0.12, 0.15, 0.65, 0.8])

# Set up the simulation object for the run we're working with, to easily
# get the relevant file paths
sim = hy.Simulation(index=sim_index)

snapshot_file = sim.get_snap_file(snap_index)
subfind_file = sim.get_subfind_file(snap_index)

# Set up a reader for the subhalo catalogue, specifying that we are only
# interested in one particular object (to speed things up). From this,
# we can access all properties of the selected subhalo as attributes
# of <subhalo>; they are loaded when first required. By default, quantities
# are returned in 'astronomically sensible' units (e.g. Mpc, Gyr, km/s).
fof_cl = hy.SplitFile(subfind_file, 'FOF', read_index=0)

# Now read in the temperature and density of the gas particles:

# Set up a reader for gas particles from the snapshot catalogue.
gas = hy.SplitFile(snapshot_file, part_type=0)

# Get all particles within 10r200 from the cluster centre
# (i.e. in the well-resolved region). Note that we could also have set up
# a ReadRegion around the cluster centre, as demonstrated in the example
# 'find_particles_near_subhalo.py' to read particles within this sphere.
ind_select = np.nonzero(
    np.linalg.norm(gas.Coordinates - fof_cl.GroupCentreOfPotential[None, :],
                   axis=1) <= 10.0*fof_cl.Group_R_Crit200)[0]

# Create the 2D histogram. The particle temperatures, densities, and
# masses are read in implicitly as they are being accessed.
histogram, xe, ye = np.histogram2d(np.log10(gas.Temperature[ind_select]),
                                   np.log10(gas.Density[ind_select]),
                                   weights=gas.Mass[ind_select],
                                   bins=[nbins, nbins],
                                   range=[temp_range, nH_range])

# Avoid NaN values in log by adding a *very* small offset to all values
histogram += sys.float_info.min

# Plot the 2D histogram
da = (xe[1]-xe[0]) * (ye[1]-ye[0])
im = plt.imshow(np.log10(histogram/da), cmap=plt.cm.inferno,
                extent=[*nH_range, *temp_range], origin='lower', aspect='auto',
                vmin=scale_range[0], vmax=scale_range[1])

ax = plt.gca()
ax.set_xlim(nH_range)
ax.set_ylim(temp_range)
ax.set_xlabel(r'$\log_{10}\,(n_\mathrm{H}\,[\mathrm{cm}^{-3}])$')
ax.set_ylabel(r'$\log_{10}\,(T\,[\mathrm{K}])$')

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
