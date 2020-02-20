"""Demonstration script to plot the phase diagram of one simulation."""

import hydrangea as hy
import numpy as np
import matplotlib.pyplot as plt

# Set script parameters
sim_index = 14                 # Which simulation do we want?
snap_index = 29                # Which snapshot do we want to analyse?
plotloc = 'phase_diagram.png'  # Where to save the output plot?

# Prepare the plot
fig = plt.figure(figsize=(5, 4))

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
# (i.e. in the well-resolved region)
ind_select = np.nonzero(
    np.linalg.norm(gas.Coordinates - fof_cl.GroupCentreOfPotential[None, :],
                   axis=1) <= 10.0*fof_cl.Group_R_Crit200)[0]

# Create the 2D histogram
histogram, xe, ye = np.histogram2d(np.log10(gas.Temperature[ind_select]),
                                   np.log10(gas.Density[ind_select]),
                                   weights=gas.Mass[ind_select],
                                   bins=[50, 50], range=[[3, 9], [-6, 3]])

da = (xe[1]-xe[0]) * (ye[1]-ye[0])

# Plot the 2D histogram
plt.imshow(np.log10(histogram/da), cmap=plt.cm.inferno,
           extent=[-6, 3, 3, 9], origin='lower')

ax = plt.gca()
ax.set_xlim((-6, 3))
ax.set_ylim((3, 9))
ax.set_xlabel(r'$\log_{10}\,(n_\mathrm{H}\,[\mathrm{cm}^{-3}])$')
ax.set_ylabel(r'$\log_{10}\,(T\,[\mathrm{K}])$')

# Save the figure
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
                        