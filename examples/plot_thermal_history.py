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
nH_range = [-9.0, 3.0]           # Plot range in (log) nH [cm^-3]
plotloc = 'thermal_history.png'  # Where to save the output plot?

# Prepare the plot
# Prepare the plot
fig = plt.figure(figsize=(5/0.8, 4))
ax1 = fig.add_axes([0.12, 0.15, 0.65, 0.8])

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

    gas = hy.SplitFile(snep_file, part_type=0)
    ind_target = np.nonzero(gas.ParticleIDs == particle_id)[0]

    # Make sure we found exactly one particle with this ID
    if len(ind_target) != 1:
        print(f"Unexpectedly found {len(ind_target)} particles with "
              f"ID {particle_id}. Something is badly wrong...")
        set_trace()

    # Assuming we found the particle exactly once, convert the
    # 1-element array to a scalar
    ind_target = ind_target[0]

    # Now that we know which index our particle is, we can set up a
    # more specific (=faster) reader for the actual properties:
    gas_part = hy.SplitFile(snep_file, part_type=0, read_index=ind_target)

    # Write out its properties into the respective arrays. They are loaded
    # in the background, only for the particle we care about. Since we have
    # not specified a unit system when setting up <gas_part>, dimensional
    # quantities are returned in 'astronomically sensible' units
    # (M_sun for mass, K for temperature, m_proton cm^-3 for densities)
    masses[isnep] = gas_part.Mass
    metallicities[isnep] = gas_part.Metallicity
    temperatures[isnep] = gas_part.Temperature
    densities[isnep] = gas_part.Densities

# Plot the quantities
plt.plot(densities, temperatures, color='grey', linewidth=1)
sc = plt.scatter(densities, temperatures, c=metallicities, vmin=0, vmax=0.02,
                 cmap=plt.cm.magma, edgecolor='lightgrey', linewidth=0.2,
                 s=masses/1.8e6*10)

ax = plt.gca()
ax.set_xlim(nH_range)
ax.set_ylim(temp_range)
ax.set_xlabel(r'$\log_{10}\,(n_\mathrm{H}\,[\mathrm{cm}^{-3}])$')
ax.set_ylabel(r'$\log_{10}\,(T\,[\mathrm{K}])$')

# Add a colour bar on the right side
ax2 = fig.add_axes([0.81, 0.15, 0.05, 0.8])
ax2.set_xticks([])
ax2.set_yticks([])
cbar = plt.colorbar(sc, cax=ax2, orientation='vertical')
fig.text(0.94, 0.55, 'Metallicity',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
