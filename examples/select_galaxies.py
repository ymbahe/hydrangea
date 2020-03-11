"""Demonstration script to select galaxies based on multiple quantities.

It illustrates how to retrieve galaxy (subhalo) properties from the
catalogues and select objects matching a set of criteria. Data are
retrieved from both the Subfind catalogue and the galaxy tables.

At the end, a simple map of the location of the selected galaxies is plotted.
"""

import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt

# Set script parameters
sim_index = 8                  # Which simulation do we want?
snap_index = 29                # Which snapshot do we want to analyse?
imsize = 10.0                  # Half-sidelength of the map, in Mpc

# Set selection criteria for subhaloes
min_subhalo_mass = 1e10        # Minimum subhalo mass to plot [M_sun]
min_stellar_mass = 5e9         # Minimum stellar mass to plot [M_sun]
max_stellar_mass = 5e10        # Maximum stellar mass to plot [M_sun]
min_m200 = 1e13                # Minimum host mass [M_sun]
satellite_status = 1           # Plot only satellites (1) or centrals (0)
min_sfr = 0.001                # Minimum SFR [M_sun/yr]
max_boundary_flag = 1          # Maximum proximity to simulation boundary

# Specify the location where to save the plot
plotloc = 'selected_subhaloes.png'

# Prepare the plot
fig = plt.figure(figsize=(4/0.8, 4))
ax1 = fig.add_axes([0.15, 0.15, 0.65, 0.8])

# Set up the simulation object for the run we're working with
sim = hy.Simulation(index=sim_index)

# To select objects matching the search criteria, we need to load the
# relevant properties for all subhaloes in the target simulation output.
# There are two ways of doing this:
# (i)  from the Subfind catalogues
# (ii) from the pre-compiled galaxy evolution tables
# Option (i) is slower, but works for all possible properties, whereas
# (ii) only works for those properties that have been pre-compiled.
# Both approaches are demonstrated here, for completeness.

# ------------ (i) Load data directly from Subfind catalogues -----------

# Form the (first) files of the required subfind and snapshot catalogues
subfind_file = sim.get_subfind_file(snap_index)

# Set up a reader for the subhalo catalogue. From this,
# we can access all properties of the subhalo as attributes
# of <subhaloes>; they are loaded when first required. By default, quantities
# are returned in 'astronomically sensible' units (e.g. Mpc, Gyr, km/s).
print("Setting up <subhaloes> reader:")
subhaloes = hy.SplitFile(subfind_file, 'Subhalo')

# Retrieve the 'boundary flag', indicating how far each subhalo is from
# the simulation edge. Values of 0 or 1 usually indicate the subhalo is
# fine, higher values should be treated with caution. This is not in the
# main Subfind catalogue, so has to be loaded from an 'extra' file.
print("Read in boundary flag:")
subhalo_extra_file = sim.sh_extra_loc
data_set = f'Snapshot_{snap_index:03d}/BoundaryFlag'
boundary_flag = hy.hdf5.read_data(subhalo_extra_file,  # Name of HDF5 file
                                  data_set)            # Name of data set

# M200 and 'satellite status' are not directly recorded in the Subfind
# catalogue, so it requires a little trickery to get them.
#
# For M200, we need to link each subhalo to its containing FOF group,
# which does have M200 recorded. To keep everyone on their toes, the
# required 'GroupNumber' property is *1-indexed*, so have to subtract 1
# to convert it to (standard) 0-indexing:
print("Find subhaloes' M200:")
fof_index = subhaloes.GroupNumber - 1
fof = hy.SplitFile(subfind_file, 'FOF')  # Set up analogous reader for FOF

# Read the M200crit values for all FOFs, and pull out the appropriate one
# for each subhalo (i.e., the one for its own FOF)
subhalo_m200 = fof.Group_M_Crit200[fof_index]

# For satellite status, the catalogue includes a property 'SubGroupNumber'
# which is 0 for centrals, and > 0 for satellites. Let's convert this to
# a binary (0 or 1) flag that can be directly compared to the
# satellite_status variable defined on top:
print("Find subhaloes' satellite status:")
subhalo_sat_flag = subhaloes.SubGroupNumber
subhalo_sat_flag[subhalo_sat_flag > 0] = 1

# Select the (indices of) subhaloes satisfying our various criteria.
# Multiple conditions can be specified by enclosing each in parentheses
# and separating them by '&' ('and', both adjoining conditions must be true
# for selection), or '|' ('or', entry is selected if at least one of the
# adjoining conditions is true). Note that np.nonzero returns a list of
# arrays, so we explicitly have to select the first (and only) entry of that
# list here (hence the final '[0]').
print("Select matching subhaloes")
sub_indices = np.nonzero((subhaloes.Mass >= min_subhalo_mass) &
                         (subhaloes.MassType[:, 4] >= min_stellar_mass) &
                         (subhaloes.MassType[:, 4] < max_stellar_mass) &
                         (subhaloes.StarFormationRate >= min_sfr) &
                         (subhalo_sat_flag == satellite_status) &
                         (subhalo_m200 >= min_m200) &
                         (boundary_flag <= max_boundary_flag))[0]

# ----------- (ii) Load data from pre-compiled galaxy catalogues ----------

# Get the catalogue file name
galaxy_catalogue = sim.fgt_loc

# Load the required properties for the target snapshot. Note that
# (i)   we use the simple HDF5 reader (hy.hdf5.read_data) here, because
#       the galaxy catalogue is stored as a single file;
# (ii)  we specify reading a particular 'column' (=snapshot) from the
#       2D arrays with the 'read_index' and 'index_dim' key words;
# (iii) most of the data is stored in log units (c.f. subhaloes above)
print("Load data from galaxy tables")
galaxy_log_mass = hy.hdf5.read_data(galaxy_catalogue, 'Msub',
                                    read_index=snap_index, index_dim=1)
galaxy_log_mstar = hy.hdf5.read_data(galaxy_catalogue, 'Mstar',
                                     read_index=snap_index, index_dim=1)
galaxy_log_sfr = hy.hdf5.read_data(galaxy_catalogue, 'SFR',
                                   read_index=snap_index, index_dim=1)
galaxy_sat_flag = hy.hdf5.read_data(galaxy_catalogue, 'SatFlag',
                                    read_index=snap_index, index_dim=1)
galaxy_log_m200 = hy.hdf5.read_data(galaxy_catalogue, 'M200',
                                    read_index=snap_index, index_dim=1)
galaxy_boundary_flag = hy.hdf5.read_data(galaxy_catalogue, 'ContFlag',
                                         read_index=snap_index, index_dim=1)

# Select the galaxy IDs satisfying our various criteria (see above for
# details on np.nonzero).
print("Select galaxies based on tables")
gal_ids = np.nonzero((galaxy_log_mass >= np.log10(min_subhalo_mass)) &
                     (galaxy_log_mstar >= np.log10(min_stellar_mass)) &
                     (galaxy_log_mstar < np.log10(max_stellar_mass)) &
                     (galaxy_log_sfr >= np.log10(min_sfr)) &
                     (galaxy_sat_flag == satellite_status) &
                     (galaxy_log_m200 >= np.log10(min_m200)) &
                     (galaxy_boundary_flag <= max_boundary_flag))[0]

# To directly compare the results to the 'subhalo-based' search above,
# we need to convert galaxy IDs to subhalo indices. For this, we need
# to select on both galaxy IDs (within the reader, using read_index
# and snap_index (outside the reader -- can only specify one inside).
# Note that we do not need to specify index_dim here, because we use
# the default (0).
print("Load galaxy to subhalo conversion index")
sub_indices_from_gal = hy.hdf5.read_data(
    galaxy_catalogue, 'SHI', read_index=gal_ids)[:, snap_index]

# Verify that the two methods give the same answer
print(f"Found {len(sub_indices)} subhaloes from Subfind, "
      f"{len(sub_indices_from_gal)} from galaxy catalogues.")

num_different = np.count_nonzero(np.sort(sub_indices)
                                 - np.sort(sub_indices_from_gal))
print(f"There are {num_different} differences between the two index arrays.")

# ------ Final part: plot the location of selected subhaloes ------------

# For simplicity, we will only do this based on the galaxy catalogue,
# using the pre-compiled positions.
position_catalogue = sim.gps_loc

# The 'Centre' data set for galaxy coordinates is 3-dimensional, and we
# need to select both on dimension 0 (galaxy IDs) and 1 (snapshot).
# As above, we specify one inside and one outside the reader.
pos_selected = hy.hdf5.read_data(position_catalogue, 'Centre',
                                 read_index=gal_ids)[:, snap_index, :]

# Convert coordinates to positions relative to the central cluster.
# To get this, we set up a special FOF reader for only the most massive
# (index-0) group. This reads only properties for that one object, which
# is much faster than re-using the (general) FOF reader set up above
fof_cl = hy.SplitFile(subfind_file, 'FOF', read_index=0)

# Note that pos_selected is a 2D array, whereas the cluster centre only 1D.
# To make numpy's indexing compare the two automatically, we must prepend
# a shallow (length-1) dimension to the latter
delta_pos = pos_selected - fof_cl.GroupCentreOfPotential[None, :]

# Plot the subhaloes (uniformly for simplicity; see the subhalo_map.py
# example for some more possibilies)
sc = plt.scatter(delta_pos[:, 0], delta_pos[:, 1],
                 marker='o', c='black', s=15)

# Some embellishments for the plot
ax1.set_xlim((-imsize, imsize))
ax1.set_ylim((-imsize, imsize))
ax1.set_xlabel(r'$\Delta x$ [Mpc]')
ax1.set_ylabel(r'$\Delta y$ [Mpc]')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False, dpi=200)
