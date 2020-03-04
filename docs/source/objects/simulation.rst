.. _simulation:

===========================================================
Automatic construction of file paths (``Simulation`` class)
===========================================================

.. |br| raw:: html

   <br />

The paths to the various simulation files can be constructed automatically, as long as the the internal file structure is mirrored from the data repository (which is automatically the case for data downloaded with the :ref:`download <download>` tool). For this, the library provides the :class:`Simulation` class.

It is instantiated by specifying the particular simulation run for which to create file paths. Simulation-wide file paths are set directly during instantiation and are available as attributes:

.. _sim_attributes:

+-----------------------+-------------------------------------------------------+
| Attribute             | Meaning                                               |
+=======================+=======================================================+
| ``run_dir``           | Base directory of the simulation                      |
+-----------------------+-------------------------------------------------------+
| ``high_level_dir``    | Directory containing high-level catalogues            |
+-----------------------+-------------------------------------------------------+
| ``fgt_loc``           | File containing pre-compiled galaxy evolution tables  |
+-----------------------+-------------------------------------------------------+
| ``gps_loc``           | File containing galaxy positions and velocities |br|  |
|                       | in all snapshots                                      |
+-----------------------+-------------------------------------------------------+
| ``galaxy_path_loc``   | File containing approximate galaxy positions and |br| |
|                       | velocities in all outputs (including snipshots)       |
+-----------------------+-------------------------------------------------------+
| ``spider_loc``        | File containing the SpiderWeb tables, which connect   |
|                       | |br| subhaloes in different snapshots to galaxies.    |
+-----------------------+-------------------------------------------------------+
| ``sh_extra_loc``      | File containing extra subhalo-level information, |br| |
|                       | beyond what is available in the Subfind catalogue     |
+-----------------------+-------------------------------------------------------+
| ``interpolation_loc`` | File containing finely interpolated orbits for |br|   |
|                       | all reasonably massive galaxies                       |
+-----------------------+-------------------------------------------------------+

For more information on these files and their content, please refer to the `data content <https://hydrangea.mpcdf.mpg.de/data_doc/index.html>`_ section of the `Hydrangea/C-EAGLE simulation website <https://hydrangea.mpcdf.mpg.de>`_. 

Examples
--------

To work with the Hydrangea version of CE-8 (including baryons, and assuming the local data repository has its root at ``/data/repository/base/``; the ``base_dir`` keyword can be ommitted if the base directory :ref:`has been set during the installation <install>`): ::

	sim = hydrangea.Simulation(8, base_dir='/data/repository/base/')

Now form the location of the galaxy tables, and the first file of the snapshot 29 particle catalogue: ::

	galaxy_tables_file = sim.fgt_loc
	snapshot_29 = sim.get_snapshot_file(29)

For the non-Hydrangea version of the same cluster: ::

	sim = hydrangea.Simulation(8, base_dir='/data/repository/base', suite_name='5r200')

For the DM-only version of CE-1 (Hydrangea), DM-only: ::
	
	sim_dm = hydrangea.Simulation(1, base_dir='/data/repository/base', sim_type='DM') 

Some more use cases are provided by the demonstration scripts in the ``examples`` directory.


Reference
---------

* :class:`hydrangea.Simulation` : Instantiate the class

Simulation-wide paths
^^^^^^^^^^^^^^^^^^^^^

* :ref:`See table above <sim_attributes>` for attributes containing simulation-wide paths.

Construct output-specific file names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* :func:`hydrangea.Simulation.get_snapshot_file` : Form a file from a snapshot catalogue
* :func:`hydrangea.Simulation.get_snipshot_file` : Form a file from a snipshot catalogue
* :func:`hydrangea.Simulation.get_subfind_file` : Form a file from a subfind catalogue

.. autoclass:: hydrangea.Simulation
   :members:
