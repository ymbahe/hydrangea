"""
Test suite to verify correct working of hydrangea utilities.
"""

import hydrangea as hy
import hydrangea.crossref as xr
import numpy as np
from pdb import set_trace

def test_crossref():
    ids_int = np.array([10, 4, 6, 2, 8, 677])
    ids_ext = np.array([8, 5, 6, 9, 2, 10, 700, 677])

    print("ids_int=", ids_int)
    print("ids_ext=", ids_ext)

    gate = xr.Gate(ids_ext, ids_int)
    ie_int, ie_matched = gate.in_int()
    
    print("ie_int=", ie_int)
    if ie_int == np.array([4, -1, 2, -1, 3, -1, 5]):
        print("Passed")
    else:
        set_trace()

    print("ie_matched=", ie_matched)
    if ie_matched == np.array([0, 2, 4, 5, 7]):
        print("Passed")
    else:
        set_trace()

    print("ids_ext[ie_matched]=", ie_ext[ie_matched])
    print("ids_int[ie_int[ie_matched]]=", ids_int[ie_int[ie_matched]])

    if ie_ext[ie_matched] == ids_int[ie_int[ie_matched]]:
        print("Passed")
    else:
        set_trace()
    
    revList = xr.create_reverse_list(ids_int)
    print("RL=", revList)
    if revList[10] == 0 and revList[4] == 1 and revList[6] == 2 and revList[2] == 3 and revList[8] == 4 and revList[677] == 5 and revList[500] == -1:
        print("Passed")
    else:
        set_trace()
    
if __name__ == "main":
    test_crossref()
    
