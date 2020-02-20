"""Demonstration script to plot the stellar mass growth of one galaxy.

This is an example of how to use information in the high-level tables
and SpiderWeb catalogue to trace the evolution of a galaxy. It also
illustrates using the read_index keyword in reading data from HDF5 files.
"""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13

# Set script parameters
sim_index = 0                # Which simulation do we want?
galaxy_id = 7679             # Which galaxy should we follow?
plotloc = 'mass_growth.png'  # Where to save the output plot?

# Prepare the plot
fig = plt.figure(figsize=(5.2, 4))

# Set up the simulation object for the run we're working with, to easily
# get the relevant file paths
sim = hy.Simulation(index=sim_index)

# ------------------------------------------------------
# I.) Simple start: the galaxy's own stellar mass growth
# ------------------------------------------------------

# Load stellar masses from the high-level catalogues (in log M/M_sun).
# By specifying a read_index, we only load the masses for the one specific
# galaxy we are interested in (at all 30 snapshots)
m_star_gal = hy.hdf5.read_data(sim.fgt_loc, 'Mstar', read_index=galaxy_id)

# The galaxy may not exist in all snapshots (either because it has not yet
# emerged above the resolution threshold, has been disrupted, or has been
# temporarily not found by Subfind). In these cases, -1 is used as a
# filler value. Let's find only snapshots with real masses:
ind_good_snap = np.nonzero(m_star_gal >= 0)[0]

# Plot the stellar mass growth history (against age of the Universe)
snap_times = hy.snep_times(time_type='age')
plt.plot(snap_times[ind_good_snap], m_star_gal[ind_good_snap],
         color='red', linewidth=2, marker='o',
         label='Main progenitor')


# ---------------------------------------------------------------------
# II.) Evolution of arbitrary quantities (not in high-level catalogues)
# ---------------------------------------------------------------------

# Now let's add a quantity that is not in the high-level catalogues
# already, such as stellar metallicity. For this, we need to go back
# to the individual subfind catalogues, so we need the subfind catalogue
# index of the galaxy in each snapshot:
shi_gal = hy.hdf5.read_data(sim.fgt_loc, 'SHI', read_index=galaxy_id)

# Now load the metallicity from each individual catalogue:
z_star = np.zeros(30)

for isnap in range(30):

    # We can skip snapshots in which the galaxy does not exist
    if shi_gal[isnap] < 0:
        continue

    # First, get the path to the (first file of the) subfind catalogue
    subfind_file = sim.get_subfind_file(isnap)

    # Set up a reader object for the multi-file catalogue (specifying
    # that we only care about one particular entry, for speed-up)
    subhalo_galaxy = hy.SplitFile(subfind_file, 'Subhalo',
                                  read_index=shi_gal[isnap],
                                  verbose=0)

    # (Stellar) metallicity is stored in a separate group 'Stars', which
    # uses a double-underscore separator to access it as an attribute
    z_star[isnap] = subhalo_galaxy.Stars__Metallicity
    
# Add the metallicity to the plot (on a second axis for clarity)
ax = plt.gca()
ax2 = ax.twinx()
plt.axes(ax2)
plt.plot(snap_times[ind_good_snap], z_star[ind_good_snap],
         color='royalblue', linewidth=2, linestyle='--')


# ---------------------------------------------------------------
# III.) Final part: total stellar mass of all progenitor galaxies
# ---------------------------------------------------------------

# So far, we have only considered the evolution of the galaxy itself
# (i.e. its main progenitor). Let's add the evolution of the total mass
# in all progenitor galaxies, including those that merge into the
# target galaxy by z = 0.

# Find all galaxies that merge with the target by z = 0: all those that
# have the target's ID as their MergeList entry
merge_ids = hy.hdf5.read_data(sim.spider_loc, 'MergeList', read_index=29,
                              index_dim=1)
ind_progenitors = np.nonzero(merge_ids == galaxy_id)[0]
print(f'Galaxy {galaxy_id} has {len(ind_progenitors)} progenitors in total.')

# We now need the stellar mass table again, but this time for all progenitor
# galaxies (in practice, we could have loaded this directly in part I)
log_m_star_progenitors = hy.hdf5.read_data(sim.fgt_loc, 'Mstar',
                                           read_index=ind_progenitors)

# Sum the progenitor masses in each snapshot (taking care of filler values)
m_star_progenitors = 10.0**log_m_star_progenitors
ind_dummy = np.nonzero(log_m_star_progenitors < 0)
m_star_progenitors[ind_dummy] = 0
m_star_tot = np.sum(m_star_progenitors, axis=0)

# Plot progenitor masses on top of the galaxy's own mass evolution
plt.axes(ax)
plt.plot(snap_times, np.log10(m_star_tot), color='indianred',
         linewidth=1, label='All progenitors')


# Some final embellishments on the figure
plt.legend(loc=4)

age_z0 = snap_times[29]
ax.set_xlim((0, age_z0))
ax.set_ylim((6, 11))
ax2.set_xlim((0, age_z0))
ax2.set_ylim((0, 0.03))

ax.set_xlabel('Age of the Universe [Gyr]')
ax.set_ylabel(r'$\log_{10}\,(M_\mathrm{star}\,/\,\mathrm{M}_\odot)$')
ax2.set_ylabel(r'Stellar metallicity ($Z$)', labelpad=20)
ax2.yaxis.label.set_rotation(-90.0)

# Colour two y axes according to which data they represent
ax2.spines['left'].set_color('red')
ax.tick_params(axis='y', colors='red')
ax.yaxis.label.set_color('red')
ax2.spines['right'].set_color('royalblue')
ax2.tick_params(axis='y', colors='royalblue')
ax2.yaxis.label.set_color('royalblue')

# Add redshifts along the upper x-axis
ax3 = ax.twiny()
zred = np.array((4,2,1,0.5,0.3,0.1,0))
age = Planck13.age(zred).value
ax3.set_xlim((0, age_z0))
ax3.set_xticks(age)
ax3.set_xticklabels([str(_z) for _z in zred])
ax3.set_xlabel('Redshift')

# Save the figure
plt.subplots_adjust(left=0.14, bottom=0.15, right=0.82, top=0.88)
plt.show
plt.savefig(plotloc, transparent=False, dpi=150)
