"""
Convenience routines for reading and writing data in HDF5 format
"""

import h5py as h5
import numpy as np
from pdb import set_trace

def check_attribute(fileName, container, attribute):
    """
    Check whether an HDF5 attribute exists.

    Parameters
    ----------
    fileName : string
        The path to the HDF5 file to test.
    container : string
        The container (group or dataset) in which the attribute's presence is
        tested (including possible containing groups).
    attribute : string
        The name of the attribute to test.

    Returns
    -------
    bool
        True if the attribute exists, False if it does not. If the container
        or file does not exist, status is also False.
    """
    try:
        f = h5.File(fileName, 'r')
        container = f[container]
        status = attribute in container.attrs.keys()
    except:
        return False
    return status

 
def check_dataset(fileName, dataset):
    """
    Check whether a specified dataset exists in an HDF5 file.
    
    Parameters
    ----------
    fileName : string
        The path to the HDF5 file to test.
    dataset : string
        The data set to test, including possible containing groups
        (e.g. 'well/nested/group/data')

    Returns
    -------
    bool
        True if fileName is an HDF5 file and contains the specified 
        data set. False otherwise.
    """
    try:
        f = h5.File(fileName, 'r')
        d = f[dataset]
    except:
        return False
    return isinstance(d, h5.Dataset)


def check_group(fileName, group):
    """
    Check whether a specified group exists in an HDF5 file.

    Parameters
    ----------
    fileName : string
        The path to the HDF5 file to test.
    group : string
        The group to test, including possible containing groups
        (e.g. 'well/nested/group/')

    Returns
    -------
    bool
        True if `fileName' is an HDF5 file and contains the specified 
        group. False otherwise.
    """
    try:
        f = h5.File(fileName, 'r')
        d = f[group]
    except:
        return False
    return isinstance(d, h5.Group)


def write_attribute(fileName, container, attName, value,
                    new=False, group=True, update=True):
    """
    Write a variable to an HDF5 attribute.

    Parameters
    ----------
    fileName : string
        The HDF5 file to write the attribute to. If it does not already exist,
        it is created new.
    container : string
        The name of the container (data set or group) to attach the attribute
        to, including possibly containing groups (e.g. '/well/nested/group/')
    attName : string
        The name of the attribute to write to.
    value : scalar, np.array, or string
        The variable to write as HDF5 attribute.
    new : bool, optional
        If a file with the specified fileName already exists, rename it to
        'fileName.old' and start a new file containing only an (empty)
        container with this attribute. Default: False.
    group : bool, optional
        If the specified container does not exist, create it as a group (True,
        default) or (empty) data set (False).
    update : bool, optional
        First check whether the attribute already exists, and update it if so 
        (default). If False, a pre-existing attribute with the same name
        is not altered, and the input value not added to the file.

    Returns
    -------
    None
    
    Note
    ----
    If the specified variable is a string, it is first converted to 
    type np.string_, in order to be acceptable as HDF5 attribute. 

    """

    # Convert None and string to HDF5-acceptable values
    if value is None:
        value = np.nan
    if isinstance(value, str):
        value = np.string_(value)

    # Move pre-existing file out of the way if desired
    if new:
        if os.path.isfile(fileName):
            os.rename(fileName, fileName+'.old')

    # Open file for creation or modification
    f = h5.File(fileName, 'a')

    # Create appropriate container if it does not yet exist
    if container not in f.keys():
        if group:
            f.create_group(container)
        else:
            f.create_dataset(container)

    # Modify existing attribute, or create new one if it does not yet exist
    if attName in f[container].attrs.keys():
        if update:
            f[container].attrs.modify(attName, value)
        else:
            print("Not changing existing attribute.")
    else:
        f[container].attrs.create(attName, value)

    # Close the HDF5 file, and we're done.
    f.close()
    return


def write_data(fileName, container, array,
               new=False, comment=None, update=True, replace=False,
               compression=None):
    """
    Write a numpy array to an HDF5 dataset.
    
    If the specified file and/or data set already exists, the default is to 
    add the data to this, modifying the pre-existing content if required.
    
    Parameters
    ----------
    fileName : string
        The HDF5 file to write the data set to. If it does not already exist,
        it is created new.
    container : string
        The name of the data set to write, including possibly containing 
        groups (e.g. '/well/nested/group/data')
    array : np.array
        The numpy array to write to an HDF5 file (can be of any data type).
    new : bool, optional
        If a file with the specified fileName already exists, rename it to
        'fileName.old' and start a new file. Default: False.
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

    Returns
    -------
    None

    """

    # Move pre-existing file out of the way if desired:
    if new:
        if os.path.isfile(fileName):
            os.rename(fileName, fileName+'.old')

    # Open file to create/modify
    f = h5.File(fileName, 'a')

    # If update is enabled, check whether data set already exists:
    if container in f:
        print("Dataset '" + container + "' already exists...")
        if not update:
            print("Not overwriting existing data.")
            return

        if replace:
            del f[container]
        else:
            dset = f[container]
            if dset.shape != array.shape:
                print("Existing dataset has wrong shape...")
                del dset
    else:
        dset = f.create_dataset(container, arr.shape, dtype = arr.dtype,
                                compression = compression)

    # Actually write array to the data set, which now definitely exists.
    dset[...] = arr

    # Write comment if desired (update it if it already exists):
    if comment is not None:
        if 'Comment' in dset.attrs.keys():
            dset.attrs.modify('Comment', np.string_(comment))    
        else:
            dset.attrs.create('Comment', np.string_(comment))

    # Close the HDF5 file, and we're done.
    f.close()
    return


def read_attribute(fileName, container, attName,
                   default=None, require=False, convert_string=True):
    """
    Read an HDF5 attribute from a specified file and container.

    Parameters
    ----------
    fileName : string
        The HDF5 file to read the attribute from. 
    container : string
        The name of the container (data set or group) containing the 
        attribute, including possibly containing groups (e.g. 'nested/group')
    attName : string
        The name of the attribute to read from the specified container.
    default : scalar, np.array, or string
        A default value to return if the attribute (or its container, or the
        file) does not exist. Default: None
    require : bool
        Raise an exception if the attribute does not exist. If False 
        (default), the 'default' value is returned instead.
    convert_string : bool
        Convert an attribute of type np.string_ to a standard string.
        Default: True

    Returns
    -------
    attribute : scalar, np.array, or string
        The attribute value read from the HDF5 file.

    """
    try:
        f = h5.File(fileName, 'r')
        att = f[container].attrs[attName]
        f.close()
    except:
        if not require:
            return default
        raise Exception("Did not find specified attribute '{:s}' in "
                        "container '{:s}'..."
                        .format(attName, container))

    if convert_string and isinstance(att, np.string_):
        att = att.decode('UTF-8')

    return att


def read_data(fileName, container, range=None, require=False):
    """
    Read one dataset from an HDF5 file.

    Optionally, only a section of the data set may be read.

    Parameters
    ----------
    fileName : string
        The HDF5 file to read the data set from. 
    container : string
        The name of the data set, including possibly containing groups 
        (e.g. 'well/nested/group/data')
    range : list or tuple of length 2, optional
        First and beyond-last element to read from the data set. If None
        (default), the entire data set is read in. If the first (last) element
        is None, the array is read from the start (to the end).
    require : bool, optional
        Require the presence of the data set, and raise an exception if it
        does not exist. If False (default), return None.

    Returns
    -------
    np.array (type and dimensions as data set).
    
    """

    try:
        f = h5.File(fileName, 'r')
        dSet = f[container]
        dType = f[container].dtype
        dShape = f[container].shape
        dSize = dShape[0]
        
        # Work out portion of data set to read
        if range is None:
            first = 0
            stop = dSize
        else:
            first = range[0]
            stop = range[1]
        if first is None:
            first = 0
        if stop is None:
            stop = dSize

        # Set up output array
        outShape = list(dShape)
        outShape[0] = stop-first
        data_out = np.empty(outShape, dType)

        # Only need to bother reading from non-zero dataset
        if dSize > 0:
            if first == 0 and stop == dSize:
                dSet.read_direct(data_out)
            else:
                data_out = dSet[first:stop, ...]
        f.close()
                
    except:
        if not require:
            return None
        raise Exception("Did not find specified data set '{:s}'..."
                        .format(container))
        
    return data_out


def attrs_to_dict(fileName, container):
    """
    Read all attributes from a specified file and container into a dict.

    Parameters
    ----------
    fileName : string
        The HDF5 file to read the attributes from. 
    container : string
        The name of the container (data set or group) containing the 
        attribute, including possibly containing groups (e.g. 'nested/group')

    Returns
    -------
    dict
        Dictionary of 'Attribute name' : 'Attribute value'

    Notes
    -----
    An exception is raised if the container does not exist. If there are
    no attributes attached to the container, an empty dict is returned.

    """

    try:
        f = h5.File(fileName, 'r')
        cont = f[container]
        attrs_tuple = list(cont.attrs.items())
    except:
        raise Exception("Could not retrieve attributes from container {:s}"
                        .format(container))
        
    att_dict = {}
    for att in attrs_tuple:
        att_dict[att[0]] = att[1]
                    
    return att_dict


def list_datasets(fileName, group=None):
    """
    Return a list of all datasets in an HDF5 group.

    Parameters
    ----------
    fileName : string
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
        f = h5.File(fileName, 'r')
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
                        .format(fileName))
    
    return ret_list


            


