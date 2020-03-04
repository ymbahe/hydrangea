.. _crossref:

=====================
Cross-matching tools
=====================

.. |br| raw:: html

   <br />

When working with simulation data, it is frequently necessary to cross-match catalogue entries in different catalogues. For example, particles in one catalogue can be matched to a different snapshot to trace the history of individual particles, or to a subfind catalogue to find particles belonging to a particular subhalo.

Although these tasks are not specific to the Hydrangea/C-EAGLE simulations, the library includes some tools for this purpose in the :mod:`hydrangea.crossref` module, both for convenience and to integrate this functionality into the :class:`SplitFile` and :class:`ReadRegion` classes.

General overview
----------------

The :mod:`hydrangea.crossref` module provides three high-level objects for cross-matching purposes:

+--------------------------------------------+----------------------------------------+
| Object                                     | Purpose                                |
+============================================+========================================+
| :class:`hydrangea.crossref.ReverseList`    | A class to compute and store a         |
|                                            | "reverse key" list, |br| which can be  |
|                                            | indexed with the key of the object     |
|                                            | |br| to be matched (e.g. a particle's  |
|                                            | ID) to find its index.                 |
+--------------------------------------------+----------------------------------------+
| :class:`hydrangea.crossref.Gate`           | A class to translate indices between   |
|                                            | two |br| catalogues according to a     |
|                                            | matched key |br| (e.g. particle IDs).  |
+--------------------------------------------+----------------------------------------+
| :func:`hydrangea.crossref.find_id_indices` | A function to find the indices         |
|                                            | corresponding to a |br| specified set  |
|                                            | of keys (e.g. particle IDs) in a       |
|                                            | |br| reference array. This is slightly |
|                                            | more direct than |br|                  |
|                                            | :class:`hydrangea.crossref.Gate` in    |
|                                            | many situations.                       |
+--------------------------------------------+----------------------------------------+

Both :class:`Gate` and :func:`find_id_indices` internally call :class:`ReverseList` when appropriate. One important limitation of the latter is that it can only be used when the keys
(typically, particle IDs) are limited to reasonably small values, because it internally sets up an array of size ``max(keys)``.

This is generally fine for all Hydrangea/C-EAGLE simulations on systems with at least 50 GB RAM, but will not work for e.g. the original EAGLE simulations
(due to a different scheme for assigning particle IDs). In such cases, :class:`Gate` and :func:`find_id_indices` switch to a different, generally slower, matching scheme based on sorting and comparing the key lists (the threshold is user-adjustable through the `max_direct` parameter in both cases).

Examples
--------

Invert a list of particle IDs (``part_ids``) to find the particle with ID 1001: ::

	rev_list = hydrangea.crossref.ReverseList(part_ids)
	index_1001 = rev_list(1001)  # Note parentheses: calling the object, not indexing an array

Match three sets of DM particle indices in snapshot 10 (``set_1``, ``set_2``, ``set_3``) to their corresponding indices in snapshot 29: ::

	# Set up particle catalogue readers at the two snapshots
	snap_10 = hydrangea.SplitFile(snapshot_file_10, 1)
	snap_29 = hydrangea.SplitFile(snapshot_file_29, 1)

	# Set up the gate
	gate = hydrangea.crossref.Gate(snap_10.ParticleIDs, snap_29.ParticleIDs)

	# Find the indices corresponding to the three sets in snapshot 29:
	ind_1, matched_1 = gate.in_int(set_1)
	ind_2, matched_2 = gate.in_int(set_2)
	ind_3, matched_3 = gate.in_int(set_3)

Here, ``ind_1``, ``ind_2``, and ``ind_3`` contain the indices in snapshot 29 corresponding to ``set_1``, ``set_2``, and ``set_3`` in snapshot 10, and ``matched_1``, ``matched_2``, ``matched_3`` the sub-indices of particles that could be matched (which should be all for DM, since all particles are present in each snapshot).

Finally, matching only a single set of particles (same set up as last example): ::

	ind_1 = hydrangea.crossref.find_id_indices(snap_10.ParticleIDs[set_1], snap_29.ParticleIDs)
	
Reference
---------

..  autoclass:: hydrangea.crossref.Gate
    :members:

..  autoclass:: hydrangea.crossref.ReverseList
    :members:

..  autofunction:: hydrangea.crossref.cKatamaran_search

..  autofunction:: hydrangea.crossref.create_reverse_list

..  autofunction:: hydrangea.crossref.find_id_indices

..  autofunction:: hydrangea.crossref.katamaran_search

..  autofunction:: hydrangea.crossref.query_array




