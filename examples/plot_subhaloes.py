"""Demonstration script to plot the location of subhaloes around a cluster."""

import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt

# Set script parameters
sim_index = 22                 # Which simulation do we want?
snap_index = 29                # Which snapshot do we want to analyse?
imsize = 25.0                  # Half-sidelength of the map, in Mpc
plotloc = 'subhalo_map.png'    # Where to save the output plot?

# Prepare the plot
fig = plt.figure(figsize=(4/0.8, 4))
ax1 = fig.add_axes([0.15, 0.15, 0.65, 0.8])

# Set up the simulation object for the run we're working with
sim = hy.Simulation(index=sim_index)

# Form the (first) files of the required subfind and snapshot catalogues
subfind_file = sim.get_subfind_file(snap_index)

# Set up a reader for the subhalo catalogue. From this,
# we can access all properties of the subhalo as attributes
# of <subhaloes>; they are loaded when first required. By default, quantities
# are returned in 'astronomically sensible' units (e.g. Mpc, Gyr, km/s).
print("Setting up <subhaloes> reader:")
subhaloes = hy.SplitFile(subfind_file, 'Subhalo')

# Also set up an equivalent reader for the FOF catalogue, to load
# various overdensity radii of the central cluster. Because we are here
# only interested in one single object, we specify a read_index to speed
# up retrieving this object's information.
fof_cl = hy.SplitFile(subfind_file, 'FOF', read_index=0)

# Extract the position of all subhaloes, relative to the cluster, in Mpc.
delta_pos = (subhaloes.CentreOfPotential
             - fof_cl.GroupCentreOfPotential[None, :])

# Retrieve the 'boundary flag', indicating how far each subhalo is from
# the simulation edge. Values of 0 or 1 usually indicate the subhalo is
# fine, higher values should be treated with caution.
subhalo_extra_file = sim.sh_extra_loc
data_set = f'Snapshot_{snap_index:03d}/BoundaryFlag'
boundary_flag = hy.hdf5.read_data(subhalo_extra_file,  # Name of HDF5 file
                                  data_set)            # Name of data set
ind_good_subhalo = np.nonzero(boundary_flag <= 1)[0]
ind_suspect_subhalo = np.nonzero(boundary_flag > 1)[0]

# Plot the subhaloes. We indicate total mass by size, and stellar mass by
# colour. Let's use different symbols for subhaloes far from and near 
# the boundary.
sc = plt.scatter(delta_pos[ind_good_subhalo, 0],
                 delta_pos[ind_good_subhalo, 1],
                 marker='o', 
                 c=np.log10(subhaloes.MassType[ind_good_subhalo, 4]),
                 s=(np.log10(subhaloes.Mass[ind_good_subhalo])-7)*7,
                 cmap=plt.cm.viridis)

plt.scatter(delta_pos[ind_suspect_subhalo, 0],
            delta_pos[ind_suspect_subhalo, 1],
            marker='X',
            c=np.log10(subhaloes.MassType[ind_suspect_subhalo, 4]),
            s=(np.log10(subhaloes.Mass[ind_suspect_subhalo])-7)*7,
            cmap=plt.cm.viridis)

ax1.set_xlim((-imsize, imsize))
ax1.set_ylim((-imsize, imsize))
ax1.set_xlabel(r'$\Delta x$ [Mpc]')
ax1.set_ylabel(r'$\Delta y$ [Mpc]')

# Add a colour bar on the right side
ax2 = fig.add_axes([0.81, 0.15, 0.05, 0.8])
ax2.set_xticks([])
ax2.set_yticks([])
cbar = plt.colorbar(sc, cax=ax2, orientation='vertical')
fig.text(0.95, 0.5, r'log$_{10}$ ($M_\mathrm{star}$ [M$_\odot$])',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False)