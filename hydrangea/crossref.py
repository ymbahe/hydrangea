import numpy as np
import ctypes as c


def katamaran_search(a, b):
    
    """
    Perform a katamaran search to locate elements of a in b.

    This assumes that the elements in a and b are unique, and
    that a and b are sorted.

    Args:
        a (array): Elements to locate. Must be unique and sorted.
        b (array): Reference array. Must also be unique and sorted.

    Returns:
        locs_a (int array): indices of a in b (-1 if not found).

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
            
        if val_a > val_b:
            ind_b += 1
            if ind_b >= len_b:
                break
            val_b = b[ind_b]
            continue
                
        if val_b > val_a:
            ind_a += 1
            if ind_a >= len_a:
                break
            val_a = a[ind_a]
            continue

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

def ckatamaran_search(a, b):
    """ 
    Sped-up version of katamaran-search which does the loop in C.

    Note that the file path to ckat.so may need updating on systems other
    than VIRGO.

    Args:
        a (array): Elements to locate. Must be unique and sorted.
        b (array): Reference array. Must also be unique and sorted.

    Returns:
        locs_a (int array): indices of a in b (-1 if not found).
    
    """

    if len(a) == 0:
        return np.zeros(0, dtype = np.int32)-1
        
    locs_a = np.zeros(len(a), dtype = np.int64)-1

    if len(b) == 0:
        return locs_a
        
    ObjectFile = "/u/ybahe/ANALYSIS/ckat.so"
    lib = c.cdll.LoadLibrary(ObjectFile)
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


def create_reverse_list(ids, delete_ids = False, cut = False, maxval = None):
    """
    Create a reverse index list from a (unique) list of IDs.

    This gives the index for each ID value, and is thus the inverse of the 
    input list (which gives the ID for an index):

    Input ID list:        --> Reverse index list:
    || 4 | 0 | 5 | 1 ||   --> || 1 | 3 | -1 | -1 | 0 | 2 ||  

    Args:
        ids (int array): Input ID list to invert. Assumed to be unique.
        delete_ids (bool, optional): Delete input ID list after building
            reverse list. Defaults to False.
        cut (bool, optional): Ignore negative input IDs. Defaults to False.
        maxval (int, optional): Build reverse list with at least maxval+1
            entries (i.e. that can be indexed by values up to maxval), 
            even if this exceeds the maximum input ID.

    Returns:
        revlist (int array): The reverse index list. 
    """

    # If there is no input, return empty array
    if len(ids) == 0 and maxval is None:
        return np.zeros(0, dtype = np.int32)
        
    if len(ids) > 0:
        maxid = ids.max()
    else:
        maxid = -1
    
    if maxval is not None:
        if maxval > maxid:
            maxid = maxval

    if len(ids) > 2e9:
        dtype = np.int64
    else:
        dtype = np.int32
    
    if cut:
        ind_good = np.nonzero(ids >= 0)[0]
        ngood = len(ind_good)
    else:
        ind_good = np.arange(ids.shape[0], dtype = int)
        ngood = ids.shape[0]

    revlist = np.zeros(np.int64(maxid+1), dtype = dtype)-1
    revlist[ids[ind_good]] = ind_good

    if delete_ids:
        del ids

    return revlist


def find_id_indices(ids, reflist, max_direct = 1e10, min_c = 1e5):
    """
    Find and return the locations of IDs in a reference list.

    The function can handle arbitrarily high ID values. If the maximum
    value is below an (adjustable) threshold, the lookup is performed
    via an explicit reverse-lookup list. Otherwise, both input and reference
    lists are sorted and compared via a "Katamaran-search".

    Note that using a reverse-lookup list is typically much faster, but 
    may use substantial amounts of memory. For small input lists, the 
    sort-and-search approach may be faster, which can be enabled by setting
    max_direct = 0.

    Args:
        ids (int array): Array of IDs whose indices should be returned.
        reflist (int array): Reference list of IDs, assumed to be unique.
            The function will search for each input ID in this array.
        max_direct (int): Maximum ID value (in input or ref list) for which
            the direct (inverse lookup list) method will be used. Defaults 
            to 1e10. Note that this will consume at least max_direct * 8 GB
            of memory.
        min_c (int): Minimum length of ID list (input or ref) for which a 
            Katamaran-search is outsourced to the C-library. Defaults to 1e5.
            Irrelevant if a direct search is performed.

    Returns:
        ind (int array): The index in reflist for each input ID (-1 if it
            could not be located at all).
        ind_match (int array): The input ID indices that could be matched.

    """
    
    maxid_in = np.max(ids)
    maxid_ref = np.max(reflist)

    if maxid_in > max_direct or maxid_ref > max_direct:
        use_direct = False
    else:
        use_direct = True

    if use_direct:
        revlist = create_reverse_list(reflist, maxval = maxid_in)
        ind = revlist[ids]
        ind_match = np.nonzero(ind >= 0)[0]

    else:
        # Need to identify matching IDs in sh_ids by brute force
        
        tstart = time.time()
        sorter_in = np.argsort(ids)
        sorter_ref = np.argsort(reflist)
        tsort = time.time()-tstart
            
        tstart = time.time()
        if len(ids) > min_c  or len(reflist) >= min_c:
            ind_in_sorted_ref = ckatamaran_search(ids[sorter_in], reflist[sorter_ref])
        else:
            ind_in_sorted_ref = katamaran_search(ids[sorter_in], reflist[sorter_ref])

        ind_prematch_in = np.nonzero(ind_in_sorted_ref >= 0)[0]
        ind = np.zeros(len(ids), dtype = int)-1
        ind[sorter_in[ind_prematch_in]] = sorter_ref[ind_in_sorted_ref[ind_prematch_in]]
        ind_match = sorter_in[ind_prematch_in]

        tkat = time.time()-tstart
        
        print("Sorting took    {:.3f} sec." .format(tsort))
        print("Kat-search took {:.3f} sec." .format(tkat))

    return ind, ind_match
