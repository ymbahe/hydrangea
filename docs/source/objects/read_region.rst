.. _read_region:

===================================================================
Reading a sub-region of a particle catalogue (``ReadRegion`` class)
===================================================================

In many cases, we only need to load properties for particles within a small, well-defined sub-volume of the simulation. :ref:`Loading the entire catalogue entry <split_file>` (with :class:`SplitFile`) is not very efficient in this case, both in terms of speed and memory usage. Instead, the :class:`ReadRegion` class can be used to only load particles within the region of interest.

..  note::
	This functionality requires that the simulation outputs are stored in a particular order, and a particle map has been set up. This is the case for all snap- and snipshots from the HYDRO runs of the Hydrangea simulations, but not for the other C-EAGLE simulations, nor for any of the DM-only counterparts. The :class:`ReadRegion` class can still be used in those cases, but it will default to reading the full catalogue entry.

General usage
-------------

The first step is to set up an instance of the class, specifying the particle type, selection volume, and possibly other parameters (see :class:`hydrangea.ReadRegion` below for all options). Once this is done, any property of the selected particles can be accessed directly as an attribute of the class instance (e.g. ``gas.Coordinates`` for the particle coordinates, assuming the instance has the name ``gas``). By default, data is returned in "astronomically sensible" units, other systems (e.g. cgs) can be specified through the ``read_units`` attribute.

For convenience, methods and properties to :ref:`read and sum particle properties directly <rr_reading>`, :ref:`connect particles to structures <rr_linking>` and :ref:`obtaining various time-stamps of the catalogue <rr_time_stamps>` are also provided.

Examples
--------

Below are a few examples of how to use the :class:`ReadRegion` class; all of them assume that ``'snapshot_file'`` :ref:`has been set up to point to one of the files <simulation>` of the snapshot catalogue to read from. More complete use cases are provided in e.g. the ``star_density.py`` and ``sf_history.py`` scripts in the ``examples`` directory.

Read all gas particle star formation rates within 100 pkpc from (100.0, 100.0, 100.0) pMpc (but accept particles at larger radii)::

	gas = hydrangea.ReadRegion(snapshot_file, 0, [100.0, 100.0, 100.0], 0.1)
	print(gas.StarFormationRate)		

As above, but using the ``exact`` keyword to only include particles strictly within the selected sphere::

	gas = hydrangea.ReadRegion(snapshot_file, 0, [100.0, 100.0, 100.0], 0.1, exact=True)

Specifying the centre *and size* of the selection region in data units (co-moving h\ :sup:`-1` Mpc) instead::

	gas = hydrangea.ReadRegion(snapshot_file, 0, [1084.0, 1084.0, 1084.0], 0.1,
	                           coordinate_units='data')

Load star particles (``part_type == 4``) within a box of size (1.0, 0.5, 3.0) pMpc, whose lower corner is at (100.0, 120.0, 110.0) pMpc::

	stars = hydrangea.ReadRegion(snapshot_file, 4, [100.0, 120.0, 110.0], [1.0, 0.5, 3.0], 
	                             shape='box', anchor_style='bottom')

.. _read_region_reference:

Reference
---------

* :class:`hydrangea.ReadRegion` : Instantiate the class

.. _rr_reading:

Reading catalogue entries
^^^^^^^^^^^^^^^^^^^^^^^^^

* :func:`hydrangea.ReadRegion.read_data` : Explicitly read a catalogue entry
* :func:`hydrangea.ReadRegion.total_in_region` : Total or average property of particles
* :func:`hydrangea.ReadRegion.get_unit_conversion` : Obtain unit conversion factors
* :attr:`hydrangea.ReadRegion.m_baryon` : Initial baryon mass
* :attr:`hydrangea.ReadRegion.m_dm` : DM particle mass

.. _rr_linking:

Linking particles to structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* :attr:`hydrangea.ReadRegion.GroupIndex` : Emulate group index catalogue entry
* :attr:`hydrangea.ReadRegion.SubhaloIndex` : Emulate subhalo index catalogue entry
* :func:`hydrangea.ReadRegion.in_subhalo` : Identify members of a particular subhalo
* :attr:`hydrangea.ReadRegion.subfind_file` : Subfind file associated to the catalogue

.. _rr_time_stamps:

Catalogue time-stamps
^^^^^^^^^^^^^^^^^^^^^

* :attr:`hydrangea.ReadRegion.aexp` : Expansion factor of catalogue
* :attr:`hydrangea.ReadRegion.lookback_time` : Lookback time in Gyr from z = 0 to the catalogue
* :attr:`hydrangea.ReadRegion.redshift` : Redshift of the catalogue
* :attr:`hydrangea.ReadRegion.time` : Age of the Universe at the catalogue


.. autoclass:: hydrangea.ReadRegion
   :members:
   :exclude-members: cantor_file, get_astro_conv, get_cgs_to_astro_factor, get_cgs_to_si_factor, get_data_to_astro_factor, get_data_to_clean_factor, get_data_to_si_factor, get_data_to_cgs_factor, get_time  
