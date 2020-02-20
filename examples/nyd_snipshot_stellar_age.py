"""Demonstration of how to use cross-matching to combine two catalogues.

The example here is how to obtain the stellar formation time in a snipshot,
for which this information is not written out. This is essentially a more
complex version of the find_particles_near_subhalo.py example, which you
may want to have a look at first.

This script also demonstrates some basic approaches for working with snipshots
in general (viz. finding the right one, finding galaxy positions)."""

# Import required packages
import numpy as np
import hydrangea as hy
import matplotlib.pyplot as plt

# Set script parameters
sim_index = 0               # Which simulation do we want?
target_z = 0.25             # Target redshift, we'll find closest snepshot
gal_index = 7679            # Which galaxy do we want?
imsize = 50                 # (Half-)size of analysis region, in kpc
nbins = 100                 # Number of bins for plotting
ptype = 4                   # Look at stars here
plotloc = 'galaxy_age.png'  # Where to save the output plot?

# First, let's find the target snepshot. Load the full list of snepshot
# times, including snipshots (including the 25-Myr-filler ones):
snep_zred = hy.snep_times(time_type='zred', snep_list='default_long')

# Find the one closest to our target redshift
ind_snep = np.argmin(np.abs(snep_zred - target_z))
print(f"Closest snepshot to target z ({target_z:.3f}) is {ind_snep}, "
      f"at z={snep_zred[ind_snep]}")

# Translate this into snepshot type and number for the particular simulation
sim = hy.Simulation(index=sim_index)
root_index, aexp, source_type, source_num = hy.get_snepshot_indices(
    sim.run_dir, 'default_long', index=ind_snep)
print("Target snepshot identified as {source_type} {source_num}.")

# Next step: find the location of the galaxy in this snepshot
gal_position = hy.hdf5.read_data(sim.galaxy_path_loc, 
                                 f'Snepshot_{root_index:04d}/Coordinates',
                                 read_index=gal_index)

# Now we can set up a ReadRegion object to load the stars close to the 
# position of the galaxy...

