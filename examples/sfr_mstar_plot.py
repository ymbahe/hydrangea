"""Example script to demonstrate reading in subhalo/galaxy properties.

Quantities are read in both from the Subfind catalogues (slow, but works for
all data) and from the pre-compiled high-level tables (faster, but only for
a subset of data).

This script demonstrates using the HDF5 utilities to read/write data
sets and attributes from/to a single file, and the SplitFile reader for
handling catalogues split over multiple files.
"""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt
import time

# Set script parameters
snap_index = 29             # Which snapshot? 29 --> z = 0
sim_index = 0               # Which simulation?
plotloc = 'sfr_mstar.png'   # Where to save the output plot?
saveloc = 'sfr_data.hdf5'   # Where to save the data?

# Set up a simulation instance to generate file paths automatically
sim = hy.Simulation(index=sim_index)


# Retrieve the 'boundary flag', indicating how far each subhalo is from
# the simulation edge. Values of 0 or 1 usually indicate the subhalo is
# fine, higher values should be treated with caution.
subhalo_extra_file = sim.sh_extra_loc
data_set = f'Snapshot_{snap_index:03d}/BoundaryFlag'
boundary_flag = hy.hdf5.read_data(subhalo_extra_file,  # Name of HDF5 file
                                  data_set)            # Name of data set

# The thresholds corresponding to each flag are stored as attributes of the
# BoundaryFlag data set, so we can retrieve them like this:
boundary_thresholds = hy.hdf5.read_attribute(
    subhalo_extra_file,      # Name of the HDF5 file
    data_set,                # Name of the data set to which attribute belongs
    'DistanceThresholds')    # Name of the attribute
print(f"Good subhaloes (code <= 1): d >= {boundary_thresholds[1]:.1f} cMpc")
print(f"Bad subhaloes (code > 3): d < {boundary_thresholds[3]:.1f} cMpc")

ind_good_subhalo = np.nonzero(boundary_flag <= 1)[0]
ind_bad_subhalo = np.nonzero(boundary_flag > 3)[0]


# Now we need to get the stellar masses and star formation rates for all
# subhaloes. The most basic way is to get them directly from the subfind
# catalogues, so let's do this first:

# Form the (first) file name of the required subfind catalogues
subfind_file = sim.get_subfind_file(snap_index)

# Set up a reader for the subhalo catalogue. From this, we can access all
# subhalo properties as attributes of <subhalo>; they are loaded when first
# required. By default, quantities are returned in 'astronomically sensible'
# units (e.g. Mpc, Gyr, km/s).
print("Setting up <subhalo> reader:")
subhalo = hy.SplitFile(subfind_file, 'Subhalo')


# Now plot the data: stellar masses are stored as the 4th column in the
# 'MassType' data set, star formation rates have their own data set. Both
# are read in when they are first accessed.
fig = plt.figure(figsize=(4, 4))

print("\n<subhalo> now reads 'MassType' and 'StarFormationRate':")
start_time = time.time()   # Time how long it takes to read the data
plt.scatter(np.log10(subhalo.MassType[ind_good_subhalo, 4]),
            subhalo.StarFormationRate[ind_good_subhalo], color='green',
            s=10, label='Clean')

print("\nNote how data are not read in again when they are re-used.")
plt.scatter(np.log10(subhalo.MassType[ind_bad_subhalo, 4]),
            subhalo.StarFormationRate[ind_bad_subhalo], color='purple',
            s=10, label='Near edge')

delta_time = time.time() - start_time
print(f"Reading the data from Subfind catalogue took {delta_time:.3f} sec.")

plt.legend()

# We can also get data in less sensible units, for example SFR in g/s. For
# this, we use the read_data function of <subhalo>. As a precaution,
# request that data be converted to 64-bit:
sfr_cgs = subhalo.read_data('StarFormationRate', units='CGS',
                            data_type=np.float64)
print(f"Subhalo 0: SFR={subhalo.StarFormationRate[0]:.2f} Msun/yr, or "
      f"{sfr_cgs[0]:.2f} g/s")

# As well as data, the redshift of the catalogue can be retrieved as
# an attribute of <subhalo>:
plt.text(11.8, 8.0, f'z = {subhalo.redshift:.2f}', va='top', ha='right')


# Put some minimal embellishments on the plot
ax = plt.gca()
ax.set_xlim((8, 12))
ax.set_ylim((1e-2, 10))
ax.set_yscale('log')

ax.set_xlabel(r'$\log_{10}\, (M_\mathrm{star}\,/\,\mathrm{M}_\odot)$')
ax.set_ylabel(r'SFR [M$_\odot$ yr$^{-1}$]')

plt.subplots_adjust(left=0.19, right=0.95, bottom=0.15, top=0.95)


# Save the figure
plt.show
plt.savefig(plotloc, transparent=False)

# The second, faster method is to use the FullGalaxyTables, which contain
# both stellar mass and SFR. Note that here the former is in log(M/M_sun),
# and SFR in log (M_sun / yr), with '-1' values for galaxies that do not
# exist in the current snapshot:

start_time = time.time()   # Again, time how long it takes to read the data

# With 'read_index' and 'index_dim', we specify that we only want to read one
# specific row in dimension 1 (data for the specified redshift).
m_star = hy.hdf5.read_data(sim.fgt_loc, 'Mstar',
                           read_index=snap_index, index_dim=1)
sfr = hy.hdf5.read_data(sim.fgt_loc, 'SFR', read_index=snap_index, index_dim=1)

delta_time = time.time() - start_time
print(f"Reading from FullGalaxyTables took {delta_time:.3f} sec.")


# Finally, let's store the extracted data in a small catalogue so we can
# do something more elaborate with it later:
hy.hdf5.write_data(saveloc, "StellarMasses", m_star,
                   comment=(f"Galaxy stellar masses from CE-{sim_index} at "
                            f"S-{snap_index}."))
hy.hdf5.write_data(saveloc, "StarFormationRates", sfr,
                   comment=f"Galaxy star formation rates.")

# Let's also write out the redshift as an attribute to StellarMasses:
hy.hdf5.write_attribute(saveloc, "StellarMasses", "Redshift",
                        subhalo.redshift)
