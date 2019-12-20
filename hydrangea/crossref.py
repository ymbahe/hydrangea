"""
Routines for matching indices between two data sets.
"""

import numpy as np
import ctypes as c
import time
import os
from pdb import set_trace

class Gate:
    """Class for linking two data sets through their particle IDs.

    It provides functionality to translate indices in one ('external') set
    into indices in another ('internal', reference) set.

    Attributes
    ----------
    revIDs : np.array(int) or None
        The reversed internal ID list (None if not created or stored).
    
    Methods
    -------
    in_int(self, index=None)
        Translates external indices to internal.
    in_ext(self, index=None)
        Translates internal indices to external.

    Note
    ----
    This class handles arbitrarily large and high input ID lists. The 
    internally stored reverse-ID list can be used to efficiently correlate
    one data set to multiple others.

    """
    
    def __init__(self, idsExt, idsInt=None, revIDs=None,
                 max_direct=int(1e10), sort_below=100, force_sort=False,
                 min_c=int(1e5), maxVal=None, save_revIDs=False,
                 verbose=False):
        """
        Parameters
        ----------
    
        idsExt : np.array(int)
            First, 'external' set of IDs.
    
        idsInt : np.array(int), optional
            Second, 'internal' set of IDs. If it is not provided (None, 
            default), the corresponding reverse-ID list must be (revIDs).
        revIDs : np.array(int), optional
            The reverse-ID list for the second, internal data set. If it is 
            not provided (None, default), the internal IDs must be (idsInt).

        max_direct : int, optional
            Maximum value in either of the ID lists for which a reverse-index 
            based match can be performed (default: 1e10). For larger IDs,
            the lists are always matched via sorting. This parameter has no 
            effect if a reverse ID list is provided.
        sort_below : int, optional
            Maximum length of either ID list for which the sort-based search 
            is always done (default: 100). Ignored if revIDs is provided.
        force_sort : bool, optional
            Force a sort-based search irrespective of input (default: False).
            This is equivalent to, but slightly faster than, setting 
            max_direct = 0. Ignored if revIDs is provided.
        min_c : int, optional
            Minimum length of longer ID list for which a (possible) 
            sort-based search is outsourced to the C-library (default: 1e5).
            This parameter has no effect if a reverse ID list is provided.
        maxVal : int or None, optional
            If a reverse list is built, make it long enough to index up 
            to value maxVal. If None (default), its length is max(idsInt)+1.
        save_revIDs : bool, optional
            Internally store the reverse ID list, if it was provided or 
            created internally (default: False). This can be useful if one 
            data set must be matched to multiple others, in which case the 
            reverse ID list can be re-used.
        verbose : bool, optional
            Provide diagnostic messages (default: False)
        """
        
        if verbose:
            print("Beginning gate setup...")

        self.revIDs = None   # Default for consistency
        
        if revIDs is None:
            # Need to set up matching instructions from scratch
            if idsInt is None:
                print("Cannot set up a gate without any info on internal "
                      "data set...")
                set_trace()

            # Can we use the direct (reverseID-list-based) method?
            use_direct = True
            if force_sort:
                use_direct = False
            elif sort_below is not None:
                if max(len(idsExt), len(idsInt)) > sort_below:
                    use_direct = False
            elif max(np.max(idsExt), np.max(idsInt)) > max_direct:
                use_direct = False

            if use_direct:
                # Quick and simple way via reverse ID list.
                # Optionally, extend reverse list for re-use purposes
                revIDs = create_reverse_list(idsInt, maxVal=maxVal)
                if save_revIDs:
                    self.revIDs = revIDs

                # Primary result is the matches of ext in int:
                self.ie_int = query_array(idsExt, revIDs)
                self.ie_matched = np.nonzero(self.ie_int >= 0)[0]
                
            else:
                # Sort-based match must be done on a case-by-case basis,
                # so can directly use find_id_indices()
                self.ie_int, self.ie_matched = find_id_indices(
                    idsExt, idsInt, force_sort=True, verbose=verbose,
                    min_c=min_c)

        else:
            # Simplest case: we were provided with a reverse list
            if save_revIDs:
                self.revIDs = revIDs

            # Account for possibility of revID array being too short
            self.ie_int = query_array(idsExt, revIDs)
            self.ie_matched = np.nonzero(self.ie_int >= 0)[0]

    def in_int(self, index=None):
        """
        Return the indices of the external data in the internal set.

        Parameters
        ----------
        index : np.array(int), optional
            Subset of the external data for which the indices should be
            returned. If None (default), indices for the whole external data 
            set are returned.

        Returns
        -------
        ie_int : np.array(int)
            The indices into the internal data set corresponding to the 
            queried external elements (-1 for unmatched external elements).
        ie_matched : np.array(int)
            The indices into the queried external elements that could be
            matched to an internal element.
        """ 

        if index is None:
            return self.ie_int, self.ie_matched
        else:
            ie_matched_index = np.nonzero(self.ie_int[index] >= 0)[0]
            return self.ie_int[index], ie_matched_index
        
    def in_ext(self, index=None):
        """
        Return the indices of the internal data in the external set.

        Parameters
        ----------
        index : np.array(int), optional
            Subset of the internal data for which the indices should be
            returned. If None (default), indices for the whole internal data 
            set are returned.

        Returns
        -------
        ii_ext : np.array(int)
            The indices into the external data set corresponding to the 
            queried internal elements (-1 for unmatched internal elements).
        ii_matched : np.array(int)
            The indices into the queried internal elements that could be
            matched to an external element.
        """ 

        ii_ext = np.zeros(self.nInt, dtype = int)-1
        ii_matched = self.ie_int[self.ie_matched]
        ii_ext[ii_matched] = np.arange(len(ii_matched))
        
        if index is None:
            return ii_ext, ii_matched
        else:
            ii_matched_index = np.nonzero(ii_ext[index] >= 0)[0]
            return ii_ext[index], ii_matched_index


def create_reverse_list(ids, delete_ids=False, assume_positive=False,
                        maxVal=None):
    """
    Create a reverse index list from a (unique) list of IDs.

    This gives the index for each ID value, and is thus the inverse of the 
    input list (which gives the ID for an index).

    Parameters
    ----------
    ids : np.array (int)
        Input ID list to invert. It is assumed to be unique, death and 
        destruction may happen if it is not.
    delete_ids : bool, optional
        Delete input ID list after building reverse list (default: False).
    assume_positive : bool, optional
        Assume that all the input values are non-negative, which speeds up
        the inversion. If False (default), the code checks explicitly 
        which input values are non-negative and only transcribes these.
    maxVal : int, optional
        Build a reverse list with at least maxval+1 entries (i.e. that can 
        be indexed by values up to maxval), even if this exceeds the maximum 
        input ID. If None (default), the reverse list is built up to 
        self-indexing capacity.

    Returns
    -------
    revList : np.array(int)
        The reverse index list. If the input list contains fewer than 
        2 billion elements, it is of type np.int32, otherwise np.int64.

    Notes
    -----
    Negative IDs in the input are ignored.

    It is a bad idea to call this function if the input ID list contains very
    large values in relation to the available memory: as a rough guide, 
    consider the limit as 1/8 * [RAM/byte]. Ignoring this will result in
    undefined, slow, and likely annoying behaviour.

    For most practical purposes, it may be more convenient to directly use
    the Gate class or the find_id_indices functions to correlate IDs, both
    of which call this function internally.

    Examples
    --------
    Standard use to invert an array:
    >>> import numpy as np
    >>> ids = np.array([4, 0, 5, 1])
    >>> create_reverse_list(ids)
    array([1, 3, -1, -1, 0, 2])

    Including negative array values:
    >>> ids = np.array([4, 0, -1, 5, 1])
    >>> create_reverse_list(ids)
    array([1, 4, -1, -1, 0, 3])

    Including use of maxVal:
    >>> create_reverse_list(ids, maxVal=10)
    array([1, 4, -1, -1, 0, 3, -1, -1, -1, -1, -1])
    """

    # If there is no input, return an empty array
    if len(ids) == 0:
        if maxVal is None:
            return np.zeros(0, dtype = np.int32)
        else:
            return np.zeros(maxVal+1, dtype = np.int32)-1

    # Find extent and required data type of output list
    maxID = ids.max()
    if maxVal is not None:
        if maxVal > maxID:
            maxID = maxVal
    if len(ids) > 2e9:
        dtype = np.int64
    else:
        dtype = np.int32

    # Do the actual inversion
    revList = np.zeros(np.int64(maxID+1), dtype = dtype)-1
    if assume_positive:
        revlist[ids] = np.arange(len(ids))
    else:
        ind_good = np.nonzero(ids >= 0)[0]
        revlist[ids[ind_good]] = ind_good

    # Delete input if necessary, to save space
    if delete_ids:
        del ids

    return revList

def find_id_indices(idsExt, idsInt, max_direct = 1e10, min_c = 1e5,
                    sort_below=None, force_sort=False, verbose=False):
    """
    Find and return the locations of IDs in a reference list.

    The function can handle arbitrarily high ID values. If the maximum
    value is below an (adjustable) threshold, the lookup is performed
    via an explicit reverse-lookup list. Otherwise, both input and reference
    lists are sorted and compared via a "Katamaran-search".

    Parameters
    ----------
    
    idsExt : np.array (int) 
        Array of IDs whose indices should be returned.
    idsInt : np.array (int)
        Internal ('reference') list of IDs, assumed to be unique. 
        The function will search for each input ID in this array.
    max_direct : int, optional
        Maximum value in either of the ID list for which a reverse-index 
        based match is performed, rather than a sort-based search 
        (default: 1e10).
    min_c : int, optional
        Minimum length of longer ID list for which a (possible) 
        Katamaran-search is outsourced to the C-library (default: 1e5).
        Ignored if a direct search is performed.
    sort_below : int or None, optional
        Maximum length of either ID list for which sort-based search is 
        preferred (default: 100). If None, always use reverseID method if 
        possible.
    force_sort : bool, optional
        Force a sorted search irrespective of maximum input values. This 
        is equivalent to, but slightly faster than, setting max_direct=0.
    verbose : bool, optional
        Print timing information; default: False

    Returns
    -------

    ie_int : np.array(int)
        The index in idsInt for each input ID (-1 if it could not be 
        located at all).
    ie_matched : np.array(int)
        The input (external) ID indices that could be matched.

    Notes
    -----
    
    For large arrays, using a reverse-lookup list is typically much faster, 
    but may use substantial amounts of memory. For small input lists, the 
    sort-and-search approach may be faster, which can be enabled with 
    force_sort=True.

    """

    # Determine whether we can use the direct (reverseID-list-based) method
    use_direct = True
    if force_sort:
        use_direct = False
    elif sort_below is not None:
        if max(len(idsExt), len(idsInt)) > sort_below:
            use_direct = False
    elif max(np.max(idsExt), np.max(idsInt)) > max_direct:
        use_direct = False

    if use_direct:
        # Quick and simple way via reverse ID list:
        revIDs = create_reverse_list(idsInt, maxVal = maxID_ext)
        ie_int = revIDs[idsExt]
        ie_matched = np.nonzero(ie_int >= 0)[0]

    else:
        # Need to identify matching IDs on an individual basis.
        tStart = time.time()
        sorterExt = np.argsort(idsExt)
        sorterInt = np.argsort(idsInt)
        if verbose:
            tSort = time.time()-tStart
            tStart = time.time()

        if max(len(idsExt), len(idsInt)) > min_c:
            # Long input list(s): use C version of Katamaran
            ind_in_sorted_int = cKatamaran_search(idsExt[sorterExt],
                                                  idsInt[sorterInt])
        else:
            # Short input lists: use python version of Katamaran
            ind_in_sorted_int = katamaran_search(idsExt[sorterExt],
                                                 idsInt[sorterInt])

            # Result is the matches of ext in int:
            subind_found = np.nonzero(ind_in_sorted_int >= 0)[0]
            ie_matched = sorterExt[subind_found]
            ie_int = np.zeros(len(idsExt), dtype = int)-1
            ie_int[ie_matched] = sorterInt[
                ind_in_sorted_int[subind_found]]

        if verbose:
            tKat = time.time()-tStart
            print("Sorting took    {:.3f} sec." .format(tSort))
            print("Kat-search took {:.3f} sec." .format(tKat))
                
    return ie_int, ie_matched

def query_array(array, indices, default=-1):
    """Retrieve values from an array that may be too short.

    Parameters
    ----------
    array : np.array
        The array to retrieve values from.
    indices : np.array(int)
        Indices of array to query.
    default : value, optional
        The value to assign to out-of-bound indices (default: -1).

    Returns
    -------
    values : np.array
        The array values of the specified indices, with out-of-bounds 
        indices set to the default.

    Note
    ----
    This function is not intended for multi-dimensional arrays, or arrays 
    of non-numerical content.

    Example
    -------
    >>> import numpy as np
    >>> arr = np.array([1, 4, 0, 8])
    >>> ind = np.array([0, 2, 8])
    >>> query_array(arr, ind, default=-100)
    array([1, 0, -100])

    """
    
    if np.max(indices) < len(array):
        return array[indices]
    else:
        values = np.zeros(len(indices), dtype = array.dtype) + default
        ind_in_range = np.nonzero(indices < len(array))[0]
        values[ind_in_range] = array[indices[ind_in_range]]
        return values
    
def katamaran_search(a, b):
    """
    Perform a "Katamaran" search to locate elements of a in b.

    This assumes that the elements in a and b are unique, and
    that a and b are sorted. 

    Parameters
    ----------
    a : np.array
        Elements to locate. Must be unique and sorted.
    b : np.array
        Reference values. Must also be unique and sorted.

    Returns
    -------
    ind_a : np.array(int)
        Indices of a in b (-1 for elements that could not be matched).

    Notes
    -----
    a and b should be of the same data type. If they contain floats, 
    the matching values must be exactly the same to machine precision.
    
    This function is implemented in pure python and may be very slow for
    large arrays. Consider using cKatamaran_search instead for these.

    """
    
    if len(a) == 0:
        return np.zeros(0, dtype = int)-1
        
    locs_a = np.zeros(len(a), dtype = int)-1
        
    if len(b) == 0:
        return locs_a
        
    ind_a = ind_b = 0
    len_a = len(a)
    len_b = len(b)
    
    val_a = a[ind_a]
    val_b = b[ind_b]

    while(True):

        # a ahead of b: test next b element
        if val_a > val_b:
            ind_b += 1
            if ind_b >= len_b:
                break
            val_b = b[ind_b]
            continue

        # b ahead of a: test next a element
        if val_b > val_a:
            ind_a += 1
            if ind_a >= len_a:
                break
            val_a = a[ind_a]
            continue

        # Match: record and increase both
        if val_a == val_b:
            locs_a[ind_a] = ind_b
            ind_a += 1
            ind_b += 1
            if ind_a >= len_a or ind_b >= len_b:
                break
            val_a = a[ind_a]
            val_b = b[ind_b]
            continue

    return locs_a

def cKatamaran_search(a, b):
    """ 
    Sped-up "Katamaran-search" which does the loop in C.

    This assumes that the elements in a and b are unique, and
    that a and b are sorted. 

    Parameters
    ----------
    a : np.array
        Elements to locate. Must be unique and sorted.
    b : np.array
        Reference values. Must also be unique and sorted.

    Returns
    -------
    ind_a : np.array(int)
        Indices of a in b (-1 for elements that could not be matched).

    Notes
    -----
    a and b should be of the same data type. If they contain floats, 
    the matching values must be exactly the same to machine precision.
    
    For small arrays, the overhead of interfacing to the C code may 
    negate the speed benefit; consider using katamaran_search here.
    
    """

    if len(a) == 0:
        return np.zeros(0, dtype = np.int32)-1
        
    locs_a = np.zeros(len(a), dtype = np.int64)-1

    if len(b) == 0:
        return locs_a

    objectDir = os.path.dirname(os.path.realpath(__file__)) + "/clib/"
    lib = c.cdll.LoadLibrary(objectDir + 'ckat.so')
    ckat = lib.ckat
    
    a_for_c = a.astype(np.int64)
    b_for_c = b.astype(np.int64)

    a_p = a_for_c.ctypes.data_as(c.c_void_p)
    b_p = b_for_c.ctypes.data_as(c.c_void_p)
        
    na_c = c.c_long(len(a))
    nb_c = c.c_long(len(b))
    locs_a_p = locs_a.ctypes.data_as(c.c_void_p)

    myargv = c.c_void_p * 5
    argv = myargv(a_p, b_p, 
                  c.addressof(na_c),
                  c.addressof(nb_c),
                  locs_a_p)
    
    succ = ckat(5, argv)
    return locs_a


