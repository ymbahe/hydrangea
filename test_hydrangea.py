"""
Test suite to verify correct working of hydrangea utilities.
"""

import hydrangea as hy
import hydrangea.crossref as xr
import numpy as np

def test_crossref():
    ids_int = np.array([10, 4, 6, 2, 8, 677])
    ids_ext = np.array([8, 5, 6, 9, 2, 10, 700, 677])

    print("ids_int=", ids_int)
    print("ids_ext=", ids_ext)

    gate = xr.Gate(ids_ext, ids_int)
    ie_int, ie_matched = gate.in_int()
    
    print("ie_int=", ie_int)
    print("ie_matched=", ie_matched)
    print("ids_ext[ie_matched]=", ie_ext[ie_matched])
    print("ids_int[ie_int[ie_matched]]=", ids_int[ie_int[ie_matched]])

    revList = xr.create_reverse_list(ids_int)
    print("RL=", revList)
    
    
if __name__ == "main":
    test_crossref()
    
