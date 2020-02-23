"""Demonstration script to plot the thermal history of a gas particle."""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt
from pdb import set_trace

# Set script parameters
sim_index = 8                    # Which simulation do we want?
particle_id = 304642637          # Which particle should we follow?
temp_range = [2.5, 9.0]          # Plot range in (log) temperature [K]
nH_range = [-8.0, 3.0]           # Plot range in (log) nH [cm^-3]
plotloc = 'thermal_history.png'  # Where to save the output plot?

# Prepare the plot
# Prepare the plot
fig = plt.figure(figsize=(5/0.8, 4))
ax1 = fig.add_axes([0.12, 0.15, 0.62, 0.8])

# Set up the simulation object for the run we're working with, to easily
# get the relevant file paths
sim = hy.Simulation(index=sim_index)

# Load the times and identities of all snepshots
root_indices, aexps, source_types, source_nums = hy.get_snepshot_indices(
    sim.run_dir, 'default_long')
n_snep = len(root_indices)
print(f"Will follow particle ID {particle_id} through {n_snep} snepshots.")

# Set up arrays to hold relevant quantities for the particle in each
# snepshot: mass, metallicity, temperature, density
masses = np.zeros(n_snep)
metallicities = np.zeros(n_snep)
temperatures = np.zeros(n_snep)
densities = np.zeros(n_snep)

# Now loop through all snepshots and find the particle's properties
for isnep in range(n_snep):

    print(f"\nProcessing snepshot {isnep}...\n")

    """
    It's worth thinking briefly about the best strategy for loading the
    properties of our particle. Reading in the full gas particle
    catalogues and picking out the right entry is rather slow, especially
    when looping through > 200 outputs.

    If we knew the index of our particle in the current catalogue, we could
    read its properties very quickly. But to find that, we would still need
    to read in the entire particle ID list.

    Instead, we can exploit the fact that our particle will not move far
    in comoving space between outputs, so we can search for it within a
    small region around its previous position. However, we don't know how
    large a search region we need, and the larger it is, the longer reading
    will take. So we start with a small region, check whether the particle
    is within it, and iteratively expand it if not. This is accomplished
    with the while loop below.

    For the first output, we have no option but to load the full ID list,
    because we do not know the location of the particle.
    """

    # Set a default search radius (in co-moving Mpc)
    search_radius = 0.1

    # Start an (indeterminate) while loop to try different search radii
    while(True):

        # Need to create the snap-/snipshot file name for current output
        if source_types[isnep] == 'snap':
            snep_file = sim.get_snapshot_file(source_nums[isnep])
        elif source_types[isnep] == 'snip':
            snep_file = sim.get_snipshot_file(source_nums[isnep])
        else:
            # If we have neither a snap- nor a snipshot, something has gone
            # badly wrong -- should not happen!
            print(f"Unexpected snepshot type '{source_types[isnep]}!")
            set_trace()

        # Can only set up a ReadRegion from output 1 onwards, for 0
        # we must load the full particle IDs
        if isnep == 0:
            gas = hy.SplitFile(snep_file, part_type=0)
        else:
            gas = hy.ReadRegion(snep_file, 0, previous_pos, search_rad,
                                coordinate_units='data')

        # (Attempt to) locate the particle in the ID list. Note that,
        # although <gas> can be either a SplitFile or ReadRegion instance,
        # data is accessed in exactly the same way for both, so we do not
        # need to distinguish between these cases here
        ind_target = np.nonzero(gas.ParticleIDs == particle_id)[0]
        
        # Check whether we found the particle
        if len(ind_target) > 0:
            break

        # We should always find the particle in the first output, because
        # there we load the full output. If we don't, it's a bad sign
        if isnep == 0:
            print("Did not find particle in first output!")
            set_trace()

        # Did not find the particle? Need to try a larger search radius
        search_rad *= 1.5

    # End of loop. Make sure we found exactly one particle with this ID        
    if len(ind_target) != 1:
        print(f"Unexpectedly found {len(ind_target)} particles with "
              f"ID {particle_id}. Something is badly wrong...")
        set_trace()

    # Assuming we found the particle exactly once, convert the
    # 1-element array to a scalar
    ind_target = ind_target[0]

    # Now that we know which index our particle is, we can set up a
    # more specific (=faster) reader for the actual properties in the first
    # output:
    if isnep == 0:
        gas_part = hy.SplitFile(snep_file, part_type=0,
                                read_index=[ind_target])
        ind_read = 0
    else:
        gas_part = gas
        ind_read = ind_target
        
    # Write out its properties into the respective arrays. They are loaded
    # in the background, only for the particle we care about. Since we have
    # not specified a unit system when setting up <gas_part>, dimensional
    # quantities are returned in 'astronomically sensible' units
    # (M_sun for mass, K for temperature, m_proton cm^-3 for densities)
    masses[isnep] = gas_part.Mass[ind_read]
    metallicities[isnep] = gas_part.Metallicity[ind_read]
    temperatures[isnep] = gas_part.Temperature[ind_read]
    densities[isnep] = gas_part.Density[ind_read]
    previous_pos = gas_part.read_data(
        'Coordinates', units='data')[ind_read, :]
    
# Plot the quantities
plt.plot(np.log10(densities), np.log10(temperatures),
         color='grey', linewidth=0.5, linestyle=':')
sc = plt.scatter(np.log10(densities), np.log10(temperatures),
                 c=metallicities, vmin=0, vmax=0.025,
                 cmap=plt.cm.magma, edgecolor='lightgrey', linewidth=0.2,
                 s=masses/1.8e6*10)
plt.scatter(np.log10(densities[-1]), np.log10(temperatures[-1]),
            facecolor='none', edgecolor='royalblue', s=30, marker='D')

ax = plt.gca()
ax.set_xlim(nH_range)
ax.set_ylim(temp_range)
ax.set_xlabel(r'$\log_{10}\,(n_\mathrm{H}\,[\mathrm{cm}^{-3}])$')
ax.set_ylabel(r'$\log_{10}\,(T\,[\mathrm{K}])$')

# Add a colour bar on the right side
ax2 = fig.add_axes([0.76, 0.15, 0.04, 0.8])
ax2.set_xticks([])
ax2.set_yticks([])
cbar = plt.colorbar(sc, cax=ax2, orientation='vertical')
fig.text(0.94, 0.55, 'Metallicity',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
set_trace()
