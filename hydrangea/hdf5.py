"""Convenience routines for reading and writing data in HDF5 format."""

import h5py as h5
import numpy as np
import os


def test_attribute(file_name, container, attribute):
    """
    Test whether an HDF5 attribute exists.

    Parameters
    ----------
    file_name : string
        The path to the HDF5 file to test.
    container : string
        The container (group or dataset) in which the attribute's presence is
        tested (including possible containing groups).
    attribute : string
        The name of the attribute to test.

    Returns
    -------
    bool
        ``True`` if the attribute exists, ``False`` if it does not. If
        the container or file does not exist, ``False`` is also returned.
    """
    try:
        f = h5.File(file_name, 'r')
        container = f[container]
        status = attribute in container.attrs.keys()
    except KeyError:
        return False
    return status


def test_dataset(file_name, dataset):
    """
    Test whether a specified dataset exists in an HDF5 file.

    Parameters
    ----------
    file_name : string
        The path to the HDF5 file to test.
    dataset : string
        The data set to test, including possible containing groups
        (e.g. 'well/nested/group/data')

    Returns
    -------
    bool
        ``True`` if file_name is an HDF5 file and contains the specified
        data set. ``False`` otherwise.
    """
    try:
        f = h5.File(file_name, 'r')
        d = f[dataset]
    except KeyError:
        return False
    return isinstance(d, h5.Dataset)


def test_group(file_name, group):
    """
    Test whether a specified group exists in an HDF5 file.

    Parameters
    ----------
    file_name : string
        The path to the HDF5 file to test.
    group : string
        The group to test, including possible containing groups
        (e.g. 'well/nested/group/')

    Returns
    -------
    bool
        ``True`` if `file_name` is an HDF5 file and contains the specified
        group. ``False`` otherwise.
    """
    try:
        f = h5.File(file_name, 'r')
        d = f[group]
    except KeyError:
        return False
    return isinstance(d, h5.Group)


def write_attribute(file_name, container, att_name, value,
                    new=False, group=True, update=True):
    """
    Write a variable to an HDF5 attribute.

    Parameters
    ----------
    file_name : string
        The HDF5 file to write the attribute to. If it does not already exist,
        it is created new.
    container : string
        The name of the container (data set or group) to attach the attribute
        to, including possibly containing groups (e.g. '/well/nested/group/')
    att_name : string
        The name of the attribute to write to.
    value : scalar, np.array, or string
        The variable to write as HDF5 attribute.
    new : bool, optional
        If a file with the specified file_name already exists, rename it to
        'file_name.old' and start a new file containing only an (empty)
        container with this attribute. Default: ``False``.
    group : bool, optional
        If the specified container does not exist, create it as a group
        (if ``True``, default) or (empty) data set (``False``).
    update : bool, optional
        First check whether the attribute already exists, and update it if so
        (default). If ``False``, a pre-existing attribute with the same name
        is not altered, and the input value not added to the file.

    Note
    ----
    If the specified variable is a string, it is first converted to
    type ``np.bytes_``, in order to be acceptable as HDF5 attribute.

    """
    # Convert None and string to HDF5-acceptable values
    if value is None:
        value = np.nan
    if isinstance(value, str):
        value = np.bytes_(value)

    # Move pre-existing file out of the way if desired
    if new:
        if os.path.isfile(file_name):
            os.rename(file_name, file_name+'.old')

    # Open file for creation or modification
    f = h5.File(file_name, 'a')

    # Create appropriate container if it does not yet exist
    if container not in f.keys():
        if group:
            f.create_group(container)
        else:
            f.create_dataset(container)

    # Modify existing attribute, or create new one if it does not yet exist
    if att_name in f[container].attrs.keys():
        if update:
            f[container].attrs.modify(att_name, value)
        else:
            print("Not changing existing attribute.")
    else:
        f[container].attrs.create(att_name, value)

    # Close the HDF5 file, and we're done.
    f.close()
    return


def write_data(file_name, container, array,
               new=False, comment=None, update=True, replace=False,
               compression=None):
    """
    Write a numpy array to an HDF5 dataset.

    If the specified file and/or data set already exists, the default is to
    add the data to this, modifying the pre-existing content if required.

    Parameters
    ----------
    file_name : string
        The HDF5 file to write the data set to. If it does not already exist,
        it is created new.
    container : string
        The name of the data set to write, including possibly containing
        groups (e.g. '/well/nested/group/data')
    array : np.array
        The numpy array to write to an HDF5 file (can be of any data type).
    new : bool, optional
        If a file with the specified file_name already exists, rename it to
        'file_name.old' and start a new file. Default: False.
    comment : string, optional
        Comment text to describe the content of the data set; written as an
        HDF5 attribute 'Comment' to the data set. If None (default), no
        comment is written.
    update : bool, optional
        Check whether the data set already exists, and update it if so
        (default: True). If False, a pre-existing data set is left intact,
        and array not written to the HDF5 file.
    replace : bool, optional
        If a data set is updated, always delete the old one first.
        If False (default), in-place update is attempted (as long as the old
        and new data sets are of the same shape). This parameter is only
        meaningful if update is True.
    compression : string, optional
        Compression to be applied to a newly created data set. Options are
        'gzip' or 'lzf' (default: None, no compression applied).
    """
    # Move pre-existing file out of the way if desired:
    if new:
        if os.path.isfile(file_name):
            os.rename(file_name, file_name+'.old')

    # Open file to create/modify
    f = h5.File(file_name, 'a')

    # If update is enabled, check whether data set already exists:
    if container in f:
        print("Dataset '" + container + "' already exists...")
        if not update:
            print("Not overwriting existing data.")
            return

        if replace:
            del f[container]

        else:
            dSet = f[container]
            if dSet.shape != array.shape:
                print("Existing dataset has wrong shape...")
                del f[container]

    # We may have just removed the data set, so have to check again
    if container not in f:
        dSet = f.create_dataset(container, array.shape, dtype=array.dtype,
                                compression=compression)

    # Actually write array to the data set, which now definitely exists.
    dSet[...] = array

    # Write comment if desired (update it if it already exists):
    if comment is not None:
        if 'Comment' in dSet.attrs.keys():
            dSet.attrs.modify('Comment', np.bytes_(comment))
        else:
            dSet.attrs.create('Comment', np.bytes_(comment))

    # Close the HDF5 file, and we're done.
    f.close()
    return


def read_attribute(file_name, container, att_name,
                   default=None, require=False, convert_string=True):
    """
    Read an HDF5 attribute from a specified file and container.

    Parameters
    ----------
    file_name : string
        The HDF5 file to read the attribute from.
    container : string
        The name of the container (data set or group) containing the
        attribute, including possibly containing groups (e.g. 'nested/group')
    att_name : string
        The name of the attribute to read from the specified container.
    default : scalar, np.array, or string
        A default value to return if the attribute (or its container, or the
        file) does not exist. Default: ``None``
    require : bool
        Raise an exception if the attribute does not exist. If False
        (default), the `default` value is returned instead.
    convert_string : bool
        Convert an attribute of type ``np.bytes_`` to a standard string.
        Default: ``True``

    Returns
    -------
    attribute : scalar, np.array, or string
        The attribute value read from the HDF5 file.

    """
    try:
        f = h5.File(file_name, 'r')
        att = f[container].attrs[att_name]
        f.close()
    except KeyError:
        if not require:
            return default
        raise Exception(f"Did not find specified attribute '{att_name}' in "
                        f"container '{container}'...")

    if convert_string and isinstance(att, np.bytes_):
        att = att.decode('UTF-8')

    return att


def read_data(file_name, container, read_range=None, read_index=None,
              index_dim=0, require=False):
    """
    Read one dataset from an HDF5 file.

    Optionally, only a section of the dataset may be read, as specified
    by the `read_range`, `read_index`, and `index_dim` parameters.

    Parameters
    ----------
    file_name : string
        The HDF5 file to read the dataset from.
    container : string
        The name of the dataset, including possibly containing groups
        (e.g. 'well/nested/group/data').
    read_range : (int, int) or ``None``, optional
        Read only elements from the first up to *but excluding* the
        second entry in the tuple (in dimension `index_dim`). If ``None``
        (default), read the entire file. Ignored if `read_index` is
        provided.
    read_index : int or np.array(int) or ``None``, optional
        Read only the specified element(s), in dimension `index_dim`. If
        `int`, a single element is read, and the first dimension truncated.
        If an array is provided, the elements between the lowest and
        highest index are read and the output then masked to the
        exact elements. If ``None`` (default), everything is read.
    index_dim : int, optional
        Dimension over which to apply the read cuts specified by
        `read_range` or `read_index` (default: 0).
    require : bool, optional
        Require the presence of the data set, and raise an exception if it
        does not exist. If ``False`` (default), return ``None`` in this
        case.

    Returns
    -------
    np.array (type and dimensions as dataset).

    """
    try:
        f = h5.File(file_name, 'r')
        dSet = f[container]
        data_type = f[container].dtype
        data_shape = f[container].shape

        # Over-write read_range if read_index is provided
        if np.isscalar(read_index):
            read_range = (read_index, read_index+1)
        elif read_index is not None:
            if len(read_index) > 0:
                read_range = (np.min(read_index), np.max(read_index) + 1)
            else:
                read_range = (0, 0)

        # Work out portion of data set to read
        out_shape = list(data_shape)
        if read_range is not None:
            first, stop = read_range
            out_shape[index_dim] = stop - first
        data_out = np.empty(out_shape, data_type)

        # Only need to bother reading non-empty data
        if f[container].size > 0 and out_shape[index_dim] > 0:
            if read_range is None:
                dSet.read_direct(data_out)
            else:
                slice_object = [np.s_[:]] * len(data_shape)
                slice_object[index_dim] = np.s_[first:stop]
                data_out = dSet[tuple(slice_object)]
        f.close()

        # Apply final truncation if needed
        if np.isscalar(read_index):
            if data_out.ndim == 1:
                data_out = data_out[0]
            else:
                slice_object = [np.s_[:]] * len(data_shape)
                slice_object[index_dim] = 0
                data_out = data_out[tuple(slice_object)]
        elif read_index is not None:
            ind_sel = read_index - read_range[0]
            slice_object = [np.s_[:]] * len(data_shape)
            slice_object[index_dim] = ind_sel
            data_out = data_out[tuple(slice_object)]

    except KeyError:
        if not require:
            return None
        raise Exception(f"Did not find specified data set '{container}'...")

    return data_out


def attrs_to_dict(file_name, container):
    """
    Read all attributes from a specified file and container into a dict.

    Parameters
    ----------
    file_name : string
        The HDF5 file to read the attributes from.
    container : string
        The name of the container (data set or group) containing the
        attribute, including possibly containing groups (e.g. 'nested/group')

    Returns
    -------
    dict
        Dictionary of 'Attribute name' : 'Attribute value'

    Note
    ----
    An exception is raised if the container does not exist. If there are
    no attributes attached to the container, an empty dict is returned.

    """
    try:
        f = h5.File(file_name, 'r')
        cont = f[container]
        attrs_tuple = list(cont.attrs.items())
    except KeyError:
        raise Exception("Could not retrieve attributes from container {:s}"
                        .format(container))

    att_dict = {}
    for att in attrs_tuple:
        att_dict[att[0]] = att[1]

    return att_dict


def list_datasets(file_name, group=None):
    """
    Return a list of all datasets in an HDF5 group.

    Parameters
    ----------
    file_name : string
        The path to the HDF5 file whose data sets should be listed.
    group : string, optional
        The group in which to look for data sets (default: None, i.e.
        search in the root group of the HDF5 file).

    Returns
    -------
    list of strings
        The names of all data sets in the target group.

    Note
    ----
    An exception is raised if the file could not be opened.
    """
    try:
        f = h5.File(file_name, 'r')
        if group is not None:
            grp = f[group]
            ret_list = [key for key in grp.keys() if isinstance(
                grp[key], h5.Dataset)]
        else:
            ret_list = [key for key in f.keys() if isinstance(
                f[key], h5.Dataset)]
        f.close()

    except:
        raise Exception("Could not read datasets from file {:s}"
                        .format(file_name))

    return ret_list
