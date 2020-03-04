.. _hdf5:

========
HDF5 I/O
========

For convenience, the library provides a number of tools to simplify reading and writing HDF5 files in the :mod:`hydrangea.hdf5` module. All of them are relatively simple wrappers around the `h5py <https://www.h5py.org/>`_ library, and can be used to read and write datasets and attributes, or inspect the contents of an HDF5 file.

..  note::
	The functions in this module are provided for coding convenience. In situations involving heavy HDF5 I/O, better performance may be obtained by using the h5py library directly.

Examples
--------

Write a numpy array to an HDF5 dataset named ``'Random'`` in a file ``data.hdf5`` (including a comment describing what the dataset contains): ::

	random_data = np.random.random(100)  # Generate some data
	hydrangea.hdf5.write_data('data.hdf5', 'Random', random_data,
	                          comment="100 randomly generated data points")

Read the data and comment back in: ::

	data_in = hydrangea.hdf5.read_data('data.hdf5', 'Random')
	data_comment = hydrangea.hdf5.read_attribute('data.hdf5', 'Comment')

Alternatively, only read data point 76: ::

	data_76 = hydrangea.hdf5.read_data('data.hdf5', 'Random', read_index=76)

Use cases more directly related to the Hydrangea/C-EAGLE simulations are provided by the demonstration scripts in the ``examples`` directory (see the :ref:`"Basic examples" <basic_examples>` section).


Reference
---------

.. _hdf5_data:

HDF5 data handling
^^^^^^^^^^^^^^^^^^

* :func:`hydrangea.hdf5.read_data` : Read an HDF5 dataset
* :func:`hydrangea.hdf5.write_data` : Write a variable to an HDF5 dataset

.. _hdf5_attributes:

HDF5 attribute handling
^^^^^^^^^^^^^^^^^^^^^^^

* :func:`hydrangea.hdf5.attrs_to_dict` : Read all HDF5 attributes from a container into a dict
* :func:`hydrangea.hdf5.read_attribute` : Read an HDF5 attribute
* :func:`hydrangea.hdf5.write_attribute` : Write a variable to an HDF5 attribute


.. _hdf5_structure:

Inspect HDF5 structure
^^^^^^^^^^^^^^^^^^^^^^

* :func:`hydrangea.hdf5.list_datasets` : List all datasets in an HDF5 group
* :func:`hydrangea.hdf5.test_attribute` : Test whether a specified HDF5 attribute exists
* :func:`hydrangea.hdf5.test_dataset` : Test whether a specified HDF5 dataset exists
* :func:`hydrangea.hdf5.test_group` : Test whether a specified HDF5 group exists


.. automodule:: hydrangea.hdf5
   :members:

