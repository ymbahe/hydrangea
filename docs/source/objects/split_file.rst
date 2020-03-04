.. _split_file:

================================================================
Read a catalogue split over multiple files (``SplitFile`` class)
================================================================

The particle and subfind catalogues of the Hydrangea/C-EAGLE simulations are distributed over several HDF5 files. Although each of these can, in principle, be read directly with `h5py <https://www.h5py.org/>`_ or the :ref:`HDF5 I/O routines <hdf5>` of this library, this is not very convenient. The :class:`SplitFile` class allows easy reading of entries from these split catalogues, without having to explicitly consider the way in which the data are distributed over files.

..  note::
	It is possible to select individual elements, or a list/range thereof, from the catalogues and read their entries much faster than the full data list. For reading particles within a particular sub-volume of a simulation, use the :ref:`ReadRegion <read_region>` class instead.

General usage
-------------

The first step is to set up an instance of the class, specifying the catalogue (e.g. particle or structure type), and possibly other parameters (see :class:`hydrangea.SplitFile` below for all options). Once this is done, any catalogue property can be accessed directly as an attribute of the class instance (e.g. ``subhalo.Mass`` for the total mass of subhaloes, assuming the instance has the name ``subhalo`` and has been set up to read a subhalo catalogue). By default, data is returned in "astronomically sensible" units, other systems (e.g. cgs) can be specified through the ``units`` attribute.

For convenience, methods and properties to :ref:`read particle properties directly <sf_reading>`, :ref:`connect particles to structures <sf_linking>` and obtaining various :ref:`time-stamps <sf_time_stamps>` and :ref:`meta-data <sf_meta>` of the catalogue.

Examples
--------

Below are a few examples of how to use the :class:`SplitFile` class to read from subfind or particle catalogues; all of them assume that ``'subfind_file'`` and ``'snapshot_file'`` :ref:`have been set up to point to one of the files <simulation>` of the subfind or snapshot catalogue, respectively, to read from. More complete use cases are provided by the demonstration scripts listed in the :ref:`"Basic examples" <basic_examples>` section.

Read the mass of all subhaloes: ::
	
	subhaloes = hydrangea.SplitFile(subfind_file, 'Subhalo')
	print(subhaloes.Mass)

Read the stellar mass (``part_type == 4``) of subhalo index 1000 (note that we do *not* specify ``part_type`` as an argument to the class constructor, because this would instruct it to read a star *particle* catalogue instead): ::

	subhalo = hydrangea.SplitFile(subfind_file, 'Subhalo', read_index=1000)
	print(subhalo.MassType[4])

Properties stored in HDF5 groups can be accessed by separating group(s) and data sets by ``'__'``. Here is how to read the stellar mass within an aperture of 30 pkpc from the subhalo centre: ::

	print(subhalo.ApertureMeasurements__Mass__030kpc[4])

If desired, entries can also be read explicitly. This also allows specifying e.g. an alternative unit system for this particular entry. Here, we read the DM mass of the subhalo in kg (why you would want to do this is another question...): ::
	
	msub_kg = subhalo.read_data('MassType[1]', units='SI', data_type='np.float64')	

Particle catalogues can be read analogously. Here, we read the formation "times" (really: expansion factors) of two specific star particles at indices 56 and 89: ::

	stars = hydrangea.SplitFile(snapshot_file, part_type=4, read_index=[56, 89])
	print(stars.StellarFormationTime)

.. _split_file_reference:

Reference
---------

* :class:`hydrangea.SplitFile` : Instantiate the class

.. _sf_reading:

Reading catalogue entries
^^^^^^^^^^^^^^^^^^^^^^^^^

* :func:`hydrangea.SplitFile.read_data` : Explicitly read a catalogue entry
* :func:`hydrangea.SplitFile.get_unit_conversion` : Obtain unit conversion factors
* :attr:`hydrangea.SplitFile.m_baryon` : Initial baryon mass
* :attr:`hydrangea.SplitFile.m_dm` : DM particle mass

.. _sf_linking:

Linking particles to structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* :attr:`hydrangea.SplitFile.GroupIndex` : Emulate group index catalogue entry
* :attr:`hydrangea.SplitFile.SubhaloIndex` : Emulate subhalo index catalogue entry
* :func:`hydrangea.SplitFile.in_subhalo` : Identify members of a particular subhalo
* :attr:`hydrangea.SplitFile.subfind_file` : Subfind file associated to the catalogue

.. _sf_time_stamps:

Catalogue time-stamps
^^^^^^^^^^^^^^^^^^^^^

* :attr:`hydrangea.SplitFile.aexp` : Expansion factor of catalogue
* :attr:`hydrangea.SplitFile.lookback_time` : Lookback time in Gyr from z = 0 to the catalogue
* :attr:`hydrangea.SplitFile.redshift` : Redshift of the catalogue
* :attr:`hydrangea.SplitFile.time` : Age of the Universe at the catalogue

.. _sf_meta:

Catalogue meta-data
^^^^^^^^^^^^^^^^^^^

* :attr:`hydrangea.SplitFile.file_offsets` : Index of first entry per file in catalogue
* :attr:`hydrangea.SplitFile.num_elem` : Number of entries to be read from catalogue
* :attr:`hydrangea.SplitFile.num_files` : Number of files in catalogue


.. autoclass:: hydrangea.SplitFile
   :members:
   :exclude-members: cantor_file, get_astro_conv, get_cgs_to_astro_factor, get_cgs_to_si_factor, get_data_to_astro_factor, get_data_to_clean_factor, get_data_to_si_factor, get_data_to_cgs_factor, get_time, read_start, read_end