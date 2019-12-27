"""
Test suite to verify correct working of hydrangea utilities.
"""

import hydrangea as hy
import hydrangea.crossref as xr
import numpy as np
from pdb import set_trace

def test_crossref():
    """Test cross-referencing modules"""

    ids_int = np.array([10, 4, 6, 2, 8, 677])
    ids_ext = np.array([8, 5, 6, 9, 2, 10, 700, 677])

    print("ids_int=", ids_int)
    print("ids_ext=", ids_ext)

    revList = xr.create_reverse_list(ids_int)
    print("RL=", revList)
    if revList[10] == 0 and revList[4] == 1 and revList[6] == 2 and revList[2] == 3 and revList[8] == 4 and revList[677] == 5 and revList[500] == -1:
        print("Passed")
    else:
        set_trace()
    
    gate = xr.Gate(ids_ext, ids_int)
    ie_int, ie_matched = gate.in_int()
    
    print("ie_int=", ie_int)
    if np.all(np.equal(ie_int, np.array([4, -1, 2, -1, 3, 0, -1, 5]))):
        print("Passed")
    else:
        set_trace()

    print("ie_matched=", ie_matched)
    if np.all(np.equal(ie_matched, np.array([0, 2, 4, 5, 7]))):
        print("Passed")
    else:
        set_trace()

    print("ids_ext[ie_matched]=", ids_ext[ie_matched])
    print("ids_int[ie_int[ie_matched]]=", ids_int[ie_int[ie_matched]])

    if np.all(np.equal(ids_ext[ie_matched], ids_int[ie_int[ie_matched]])):
        print("Passed")
    else:
        set_trace()

    ii_ext, ii_matched = gate.in_ext()
    print("ii_ext=", ii_ext)
    if np.all(np.equal(ii_ext, np.array([5, -1, 2, 4, 0, 7]))):
        print("Passed")
    else:
        set_trace()

    print("ii_matched=", ii_matched)
    

    inds_fii,  matched_fii = xr.find_id_indices(ids_ext, ids_int, min_c = 0)
    print("inds_fii=", inds_fii)
    if np.all(np.equal(inds_fii, np.array([4, -1, 2, -1, 3, 0, -1, 5]))):
        print("Passed")
    else:
        set_trace()

    print("matched_fii=", matched_fii)
    if np.all(np.equal(matched_fii, np.array([0, 2, 4, 5, 7]))):
        print("Passed")
    else:
        set_trace()


    subind_ext = np.array([3, 1, 0])
    subind_int, submatched = gate.in_int(subind_ext)
    print("subind_ext=", subind_ext)
    print("subind_int=", subind_int)
    print("submatched=", submatched)
    if np.all(np.equal(subind_int, np.array([-1, -1, 4]))):
        print("Passed")
    else:
        set_trace()
    if np.all(np.equal(submatched, np.array([2]))):
        print("Passed")
    else:
        set_trace()

    subind_int = np.array([4, 1, 5, 5])
    subind_ext, submatched = gate.in_ext(subind_int)
    print("subind_int=", subind_int)
    print("subind_ext=", subind_ext)
    print("submatched=", submatched)
    if np.all(np.equal(subind_ext, np.array([0, -1, 7, 7]))):
        print("Passed")
    else:
        set_trace()
    if np.all(np.equal(submatched, np.array([0, 2, 3]))):
        print("Passed")
    else:
        set_trace()

    test = xr.query_array(ids_ext, np.array([0, 4, 3, 2, -1, 100]),
                          default = -10)
    print("test=", test)
    if np.all(np.equal(test, np.array([8, 2, 9, 6, -10, -10]))):
        print("Passed")
    else:
        set_trace()


def check_equal(arr1, arr2):
    if np.all(np.equal(arr1, arr2)):
        print("Passed")
    else:
        set_trace()
        
if __name__ == "__main__":
    test_crossref()
    
