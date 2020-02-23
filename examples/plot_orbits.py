"""Demonstration script to plot the orbit of one galaxy in HYDRO and DM.

It demonstrates reading data from the pre-compiled galaxy tables, matching
subhaloes between the HYDRO and DM-ONLY versions of a simulation, and 
working with the interpolated high-time-resolution orbit tables.
"""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt
import matplotlib.colors as colors   # For truncated grey-scale color map
from pdb import set_trace

# Set script parameters
sim_index = 0                  # Which simulation do we want?
galaxy_id = 123                # Galaxy ID (in HYDRO) to plot
plotloc = 'orbit_plot.png'     # Where to save the output plot?

# Set up the simulation object for the run we're working with, to easily
# get the relevant file paths
sim_hydro = hy.Simulation(index=sim_index)
sim_dm_only = hy.Simulation(index=sim_index, sim_type='DM')

# Prepare the plot
fig = plt.figure(figsize=(4*0.8/0.6, 4))
ax1 = fig.add_axes((0.15, 0.15, 0.65, 0.8))

# ------------------------------------------------------------------
# Part I: Find the corresponding galaxy ID in the DM-only simulation
# ------------------------------------------------------------------

# Find the corresponding galaxy ID in the DM simulation. For this, we find
# the snapshot in which the galaxy has its maximum total mass, and then
# match its subhalo
mtot_hydro = hy.hdf5.read_data(sim_hydro.fgt_loc, 'Msub',
                               read_index=galaxy_id)
snap_max = np.argmax(mtot_hydro)
print(f"Galaxy {galaxy_id} reaches its max subhalo mass in snap {snap_max}.")

# Now we find the subhalo index of the galaxy in this snapshot
shi_hydro_max = hy.hdf5.read_data(
    sim_hydro.fgt_loc, 'SHI', read_index=galaxy_id)[snap_max]

# And now we find the matching DM-only subhalo in the same snapshot
shi_dm_max = hy.hdf5.read_data(sim_hydro.sh_extra_loc,
                               f'Snapshot_{snap_max:03d}/MatchInDM',
                               read_index=shi_hydro_max)

# And find the DM-only galaxy ID
galaxy_id_dmo = hy.hdf5.read_data(sim_dm_only.spider_loc,
                                  f'Subhalo/Snapshot_{snap_max:03d}/Galaxy',
                                  read_index=shi_dm_max)
print(f"Galaxy {galaxy_id} in HYDRO matched to {galaxy_id_dmo} in DM-only")


# ---------------------------------------------------------------
# Part II: Extract the orbits of the galaxies in both simulations
# ---------------------------------------------------------------

# Find the central of the galaxy at z = 0 in both simulations
central_id_hydro = hy.hdf5.read_data(sim_hydro.fgt_loc, 'CenGal',
                                     read_index=galaxy_id)[29]
central_id_dmo = hy.hdf5.read_data(sim_dm_only.fgt_loc, 'CenGal',
                                   read_index=galaxy_id_dmo)[29]


# We now extract the high-time-resolution orbits of all four galaxies
def get_galaxy_paths(galaxy_id, sim, load_times=False):
    """Extract the galaxy paths for a specified galaxy in a simulation."""
    # First, we must look up the galaxy in the interpolation list
    interp_id = hy.hdf5.read_data(sim.interpolation_loc, 'GalaxyRevIndex',
                                  read_index=galaxy_id)

    # If the galaxy is not in the interpolation list, quit
    if interp_id < 0:
        print(f"Could not locate galaxy {galaxy_id} in interpolation list.")
        return

    # Extract the interpolated orbit
    orbit_pos = hy.hdf5.read_data(sim.interpolation_loc,
                                  'InterpolatedPositions',
                                  read_index=interp_id)

    # If needed, also load the interpolation times
    if not load_times:
        return orbit_pos

    orbit_times = hy.hdf5.read_data(sim.interpolation_loc,
                                    'InterpolationTimes')
    return orbit_pos, orbit_times


orbit_hydro, orbit_times = get_galaxy_paths(galaxy_id, sim_hydro,
                                            load_times=True)
orbit_hydro -= get_galaxy_paths(central_id_hydro, sim_hydro)
orbit_dmo = (get_galaxy_paths(galaxy_id_dmo, sim_dm_only)
             - get_galaxy_paths(central_id_dmo, sim_dm_only))


# -------------------------
# Part III: Plot the orbits
# -------------------------

# Create a truncated grey color map (Greys_r becomes white towards 1)
# [see https://stackoverflow.com/questions/18926031/]
cmap_greys = colors.LinearSegmentedColormap.from_list(
    'trunc(Greys_r, 0.0, 0.8)', plt.cm.Greys_r(np.linspace(0.0, 0.8, 1000)))

plt.sca(ax1)
sc_dm = plt.scatter(orbit_dmo[0, :], orbit_dmo[1, :], c=orbit_times,
                    edgecolor='none', vmin=0, vmax=14.0, cmap=cmap_greys,
                    s=10)
sc_hy = plt.scatter(orbit_hydro[0, :], orbit_hydro[1, :], c=orbit_times,
                    edgecolor='none', vmin=0, vmax=14.0, cmap=plt.cm.viridis,
                    s=10)

# Some embellishments on the plot
ax1.set_xlim((-0.8, 0.8))
ax1.set_ylim((-0.8, 0.8))
ax1.set_xlabel(r'$\Delta x$ [pMpc]')
ax1.set_ylabel(r'$\Delta y$ [pMpc]')

# Add colour bars on the side
ax2 = fig.add_axes([0.81, 0.15, 0.03, 0.8])
ax2.set_xticks([])
ax2.set_yticks([])
cbar_dm = plt.colorbar(sc_dm, cax=ax2, orientation='vertical')
cbar_dm.set_ticks([])

ax3 = fig.add_axes([0.84, 0.15, 0.03, 0.8])
ax3.set_xticks([])
ax3.set_yticks([])
cbar_hy = plt.colorbar(sc_hy, cax=ax3, orientation='vertical')

fig.text(0.94, 0.55, 'Age of the Universe [Gyr]',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
