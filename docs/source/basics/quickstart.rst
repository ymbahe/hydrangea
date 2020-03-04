.. _five_key_elements:

=================
Five key elements
=================

This section lists and briefly describes five of the most important elements of the hydrangea library. For a full description of all elements and their options, including 
those mentioned below, please see the detailed documentation pages below.

Download data
=============

Coming soon...


Auto-generate file names
========================

The :class:`hydrangea.Simulation` class can be used to automatically generate the file names of the various catalogues. For example, to get the high-level galaxy catalogue file, and the first file of the snapshot 26 (`z = 0.101`) particle and subfind catalogues for simulation CE-29: ::

	sim = hydrangea.Simulation(29)
	galaxy_file = sim.fgt_loc
	snapshot_file = sim.get_snapshot_file(26)
	subfind_file = sim.get_subfind_file(26)

For more information and details, :ref:`see here <simulation>`.

Read single HDF5 files
======================

The :func:`hdf5.read_data` function provides a convenient wrapper around the `h5py <https://www.h5py.org/>`_ library, to read data from an HDF5 file to a numpy array in a single line. For example, to read the stellar masses of all galaxies in all snapshots of simulation CE-29 (in `log`\ :sub:`10` [`M`/M\ :sub:`Sun`\ ]): ::

	galaxy_log_mstar = hydrangea.hdf5.read_data(galaxy_file, 'Mstar')

For more information and details, :ref:`see here <hdf5>`.


Read catalogues
================

The :class:`hydrangea.SplitFile` class can be used to read in full snapshot (particle) catalogues, or substructure catalogues from Subfind. In both cases, it is also possible to read only entries at a given index, or a list (or range) of indices. For example, to read the (total) initial stellar masses of all subhaloes in snapshot 26 (in M\ :sub:`Sun`): ::

	subhaloes = hydrangea.SplitFile(subfind_file, 'Subhalo')
	print(subhaloes.StellarInitialMass)  # In practice, would do something more meaningful

More information and details are given :ref:`here <split_file>`.


Read sub-regions of snapshots
=============================

With the :class:`ReadRegion` class, a specified sub-region of a snapshot can be read much more quickly than the full particle catalogue. For example, to read the temperatures of all gas particles within 100 pkpc (= 0.1 pMpc) from the centre of potential of subhalo 0 in snapshot 26 (in K): ::

	subhalo_0 = hydrangea.SplitFile(subfind_file, 'Subhalo', read_index=0)
	gas = hydrangea.ReadRegion(snapshot_file, 0, subhalo_0.Coordinates, 0.1, exact=True)
	print(gas.Temperature)  # In practice, would do something more meaningful

More information and details are given :ref:`here <read_region>`.

