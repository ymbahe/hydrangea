.. _examples:

===============
Example scripts
===============

.. |br| raw:: html

   <br />

The ``examples`` directory contains a number of example python scripts that demonstrate ways of using (elements of) this library for various science analysis steps. Each script is a self-contained unit which produces a plot from the simulation data. The scripts are heavily commented to highlight their internal workings. It is a good idea to run and modify (some of) them yourself before using the library for your own scientific goals!

..  note::
		
	The emphasis of these examples is to demonstrate how to use various parts of the library. The science questions addressed by them are selected to fit with the library aspects illustrated, rather than representing key scientific targets of the simulations. Also note that the plotting part is deliberately kept simple, so these examples should not be regarded as recommendations on how to present scientific data.
	
.. _basic_examples:

Basic examples
--------------

These are relatively straightforward, and probably the best ones to try out first.

+---------------------------+-------------------------------------------+--------------------------+
| Name                      | Library/simulation features covered       | Science aspect           |
+===========================+===========================================+==========================+
| :ref:`x-sfr-mstar-plot`   | - Read/write data from/to HDF5 files      | Plot M\ :sub:`star` -    |
|                           | - Read (full) subhalo catalogue |br|      | SFR relation             | 
|                           |   (:class:`SplitFile`)                    |                          |
|                           | - Avoiding the simulation boundary        |                          |
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-subhalo-map`      | Use :class:`SplitFile` to                 | Plot subhaloes |br|      |
|                           |                                           | around a cluster         |
|                           | - read an individual FOF catalogue entry  |                          |
|                           | - read a (full) subhalo catalogue         |                          | 
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-star-density`     | - Read data for one subhalo |br|          | Plot galaxy image        |
|                           |   (:class:`SplitFile`)                    |                          |
|                           | - Read particles in a small volume |br|   |                          |
|                           |   (:class:`ReadRegion`)                   |                          |
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-phase-diagram`    | - Use :class:`~hydrangea.SplitFile` to    | Plot the gas phase |br|  |
|                           |                                           | diagram of a simulation  |
|                           |   - read an entire particle catalogue     |                          |
|                           |   - read an individual entry from the     |                          |
|                           |     |br| FOF catalogue                    |                          |
|                           |                                           |                          |
|                           | - Avoiding the simulation boundary        |                          |
+---------------------------+-------------------------------------------+--------------------------+

..  toctree::
    :hidden:

    examples/sfr_mstar_plot
    examples/subhalo_map
    examples/star_density
    examples/phase_diagram


.. _complex_examples:


More complex example scripts
-----------------------------

These examples are slightly longer and/or more complex than the "basic" examples. In the second column, only features that have not already been covered in the basic examples above are listed.

+---------------------------+-------------------------------------------+--------------------------+
| Name                      | Library/simulation features covered       | Science aspect           |
+===========================+===========================================+==========================+
| :ref:`x-galaxy-stars`     | :func:`in_subhalo()` function to |br|     | Plot images of a |br|    |
|                           | test particles for subhalo membership     | disrupting galaxy        | 
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-galaxy-evolution` | - Follow one galaxy's evolution           | Plot M\ :sub:`star` and  |
|                           | - Find all progenitors of a galaxy |br|   | Z\ :sub:`star` |br|      | 
|                           |   in a given snapshot                     | evolution for a galaxy   |
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-sf-history`       | - Use :class:`ReadRegion` to load         | Plot SF history of |br|  |
|                           |   particles |br| within the resolved      | central and satellite    |
|                           |   simulation region                       | |br| galaxies            |
|                           | - Virtual ``SubhaloIndex``                |                          | 
|                           |   element to |br| connect particles to    |                          |
|                           |   their subhaloes                         |                          |
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-snipshot-age`     | - Working with snipshots                  | Plot stellar age         |
|                           | - :class:`Gate` class for matching        | map |br| of a galaxy     |
|                           |   entries |br| between two catalogues     | at a |br| specified time |
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-trace-particles`  | :func:`find_id_indices` function to       | Phase diagram of gas     |
|                           | match |br| entries between catalogues     | |br| destined to sit in  |
|                           | |br| (alternative to :class:`Gate`        | the |br| cluster core    |
|                           | class)                                    |                          |
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-thermal-history`  | - Use full set of snip- and snapshots     | Path of one particle     |
|                           | - Alternative units in :class:`SplitFile` | |br| across the phase    |
|                           |   |br| and :class:`ReadRegion` readers    | |br| diagram over time   |
+---------------------------+-------------------------------------------+--------------------------+
| :ref:`x-orbits`           | - Match subhaloes between HYDRO and |br|  | Compare the orbit of     |
|                           |   DM-only simulations                     | |br| a galaxy in HYDRO   |
|                           | - Interpolated orbits for very high |br|  | |br| and DM-only runs    |
|                           |   time resolution                         |                          |
+---------------------------+-------------------------------------------+--------------------------+

..  toctree::
    :hidden:

    examples/galaxy_stars
    examples/galaxy_evolution
    examples/sf_history
    examples/snipshot_age
    examples/trace_particles
    examples/thermal_history
    examples/orbits
