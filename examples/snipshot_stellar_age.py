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
from pdb import set_trace

# Set script parameters
sim_index = 0               # Which simulation do we want?
target_z = 0.25             # Target redshift, we'll find closest snepshot
gal_id = 7679               # Which galaxy do we want?
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
                                 read_index=gal_id)

# Now we can set up a ReadRegion object to load the stars close to the 
# position of the galaxy...
print("\nSetting up <parts> reader:")
if source_type == 'snap':
    snepshot_file = sim.get_snap_file(source_num)
elif source_type == 'snip':
    snepshot_file = sim.get_snip_file(source_num)
else:
    print(f"Unexpected source type '{source_type}'!")
    set_trace()

parts = hy.ReadRegion(snepshot_file,      # (First) file to read
                      ptype,              # Particle type to read
                      gal_position,       # Centre of region
                      imsize/1e3,         # Size of region [pMpc]
                      shape='cube')       # Shape of region
print("")

# To get the stellar age, we need to link each particle to a subsequent
# snapshot.
ref_snap_index = hy.get_next_snapshot(parts.aexp)
ref_snap_file = sim.get_snap_file(ref_snap_index)
print("Next snapshot is {ref_snap_index}.")

# For simplicity, we will only load particles around the target galaxy.
# In principle, we could also set up a SplitFile reader to load the entire
# star catalogue in the reference snapshot.
gal_position_ref = hy.hdf5.read_data(
    sim.gps_loc, 'Centre', read_index=gal_id)[ref_snap_index, :]
parts_ref = hy.ReadRegion(ref_snap_file, ptype, gal_position_ref, imsize/1e2)

# Now set up a 'gate' to cross-match particles between the target snepshot
# and the reference snapshot:
gate = hy.crossref.Gate(parts.ParticleIDs, parts_ref.ParticleIDs)

# We can now access the indices of particles in the (selected region of the)
# reference snapshot corresponding to indices in the target snipshot:
indices_in_ref_all = gate.in_int()
snip_formation_aexp = parts_ref.StellarFormationTime[indices_in_ref_all]

# Aside: if we only wanted the (ref-snapshot) indices for e.g. the first 100
# particles in the target snepshot read-region, we could get these directly
# by specifying them in the call to in_int():
indices_in_ref_first = gate.in_int(np.arange(100))

# Translate formation aexp to age at the time of the target snepshot
# (aexp_to_time gives the lookback time to z=0, from which we subtract the
# cosmic time at the target snepshot)
snip_stars_age = (hy.aexp_to_time(snip_formation_aexp, 'lbt')
                  - parts.time)

# Get particle coordinates relative to subhalo centre (in kpc)
print("\nNote how <parts> now reads the 'Coordinates' and 'Mass' data sets:")
part_relpos = (parts.Coordinates - gal_position[None, :]) * 1e3

# Make a simple 2D histogram map of mass-weighted stellar age
# Note that we here don't convert the (temporary) mass_map to proper
# surface density units (c.f. find_particles_near_subhalo.py example),
# because we only need it to normalise the age map
mass_map = np.histogram2d(part_relpos[:, 1], part_relpos[:, 0], bins=nbins,
                          range=[[-imsize, imsize], [-imsize, imsize]],
                          weights=parts.Mass)[0]
age_map = np.histogram2d(part_relpos[:, 1], part_relpos[:, 0], bins=nbins,
                         range=[[-imsize, imsize], [-imsize, imsize]],
                         weights=parts.Mass*snip_stars_age)[0] / mass_map

# Plot the age map, in log scale
fig = plt.figure(figsize=(4/0.8, 4))
ax1 = fig.add_axes([0.15, 0.15, 0.65, 0.8])
log_age = np.log10(age_map)
im = plt.imshow(log_age,
                extent=[-imsize, imsize, -imsize, imsize],
                aspect='equal', interpolation='nearest',
                origin='lower', alpha=1.0, cmap=plt.cm.inferno,
                vmin=-1, vmax=10)

ax1.set_xlim((-imsize, imsize))
ax1.set_ylim((-imsize, imsize))
ax1.set_xlabel(r'$\Delta x$ [kpc]')
ax1.set_ylabel(r'$\Delta y$ [kpc]')

# Add a colour bar on the right side
ax2 = fig.add_axes([0.81, 0.15, 0.05, 0.8])
ax2.set_xticks([])
ax2.set_yticks([])
cbar = plt.colorbar(im, cax=ax2, orientation='vertical')
fig.text(0.95, 0.5, r'log$_{10}$ (Age [Gyr])',
         rotation=90.0, va='center', ha='left')

# Save the figure
plt.show
plt.savefig(plotloc, transparent=False)
