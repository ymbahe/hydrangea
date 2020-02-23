"""Demonstration script to plot the orbit of one galaxy in HYDRO and DM."""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt

# Set script parameters
sim_index = 0                  # Which simulation do we want?
galaxy_id = 123                # Galaxy ID (in HYDRO) to plot
plotloc = 'orbit_plot.png'     # Where to save the output plot?

# Set up the simulation object for the run we're working with, to easily
# get the relevant file paths
sim_hydro = hy.Simulation(index=sim_index)
sim_dm_only = hy.Simulation(index=sim_index, sim_type='DM')

# Prepare the plot
fig = plt.figure(figsize=(4, 4))


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
    interp_id = hy.hdf5.read_data(sim.interp_loc, 'GalaxyRevIndex',
                                  read_index=galaxy_id)

    # If the galaxy is not in the interpolation list, quit
    if interp_id < 0:
        print(f"Could not locate galaxy {galaxy_id} in interpolation list.")
        return

    # Extract the interpolated orbit
    orbit_pos = hy.hdf5.read_data(sim.interp_loc, 'InterpolatedPositions',
                                  read_index=interp_id)

    # If needed, also load the interpolation times
    if not load_times:
        return orbit_pos

    orbit_times = hy.hdf5.read_data(sim.interp_loc, 'InterpolationTimes')
    return orbit_pos, orbit_times


orbit_hydro, orbit_times = get_galaxy_paths(galaxy_id, sim_hydro,
                                            load_times=True)
orbit_hydro -= get_galaxy_paths(central_id_hydro, sim_hydro)
orbit_dmo = (get_galaxy_paths(galaxy_id_dmo, sim_dm_only)
             - get_galaxy_paths(central_id_dmo, sim_dm_only))


# -------------------------
# Part III: Plot the orbits
# -------------------------

plt.plot(orbit_hydro[0, :], orbit_hydro[1, :], color='green', linewidth=0.5)
plt.scatter(orbit_hydro[0, :], orbit_hydro[1, :], c=orbit_times,
            edgecolor='none', vmin=0, vmax=14.0, cmap=plt.cm.viridis)

plt.plot(orbit_dmo[0, :], orbit_dmo[1, :], color='black', linewidth=0.5)
plt.scatter(orbit_dmo[0, :], orbit_dmo[1, :], c=orbit_times,
            edgecolor='none', vmin=0, vmax=14.0, cmap=plt.cm.Greys)


# Some embellishments on the plot
ax = plt.gca()
ax.set_xlim((-1, 1))
ax.set_ylim((-1, 1))
ax.set_xlabel(r'$\Delta x$ [Mpc]')
ax.set_ylabel(r'$\Delta y$ [Mpc]')

# Save the figure
plt.subplots_adjust(left=0.14, bottom=0.15, right=0.82, top=0.88)
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
