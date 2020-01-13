"""
Test suite to verify correct working of hydrangea utilities.
"""

import hydrangea as hy
import hydrangea.hdf5 as hd
import hydrangea.crossref as xr
import numpy as np
from pdb import set_trace
import sim_tools as st

rundir = '/virgo/simulations/Hydrangea/10r200/CE-11/HYDRO/'
isnap = 22

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

def test_form_files():
    snapdir, subdir = hy.form_files(rundir, isnap, 'snap sub')
    magdir = hy.form_files(rundir, isnap, 'partmag')
    print("Please check following names:")
    print(snapdir, subdir, magdir)
    set_trace()
    
def test_read_region():
    snapdir = hy.form_files(rundir, isnap, 'snap')
    posloc = rundir + '/highlev/GalaxyPositionsSnap.hdf5'
    galpos = hd.read_data(posloc, 'Centre', range = (299, 300))[0, isnap, :]
    if galpos is None:
        set_trace()
    
    gpos_full = st.eagleread(snapdir, 'PartType0/Coordinates', astro = True)[0]
    gids_full = st.eagleread(snapdir, 'PartType0/ParticleIDs', astro = False)
    gmass_full = st.eagleread(snapdir, 'PartType0/Mass', astro = True)[0]

    # TEST SPHERE #
    readReg = hy.ReadRegion(snapdir, 0, [*galpos-1.0, 2.0], shape='Sphere',
                            verbose=True, exact=True)
    gpos = readReg.read_data('Coordinates')
    gids = readReg.read_data('ParticleIDs')
    grad = np.linalg.norm(gpos-galpos[None, :], axis = 1)
    ind_sphere = np.nonzero(grad <= 1.0)[0]

    if len(ind_sphere) == len(grad):
        print("Passed")
    else:
        set_trace()

    mtot = readReg.total_in_region('Mass')

    grad_full = np.linalg.norm(gpos_full - galpos[None, :], axis = 1)
    ind_sphere_full = np.nonzero(grad_full <= 1.0)[0]
    if len(ind_sphere_full) == len(ind_sphere):
        print("Passed")
    else:
        set_trace()

    if np.sum(gmass_full[ind_sphere_full]) == mtot:
        print("Passed")
    else:
        set_trace()

    gv = readReg.total_in_region('Velocity', weightQuant='Mass', astro=True)

    gvel = readReg.read_data('Velocity', astro=True)
    gmass = readReg.read_data('Mass', astro=True)
    gv_check = np.sum(gvel*gmass[:, np.newaxis], axis = 0)/np.sum(gmass)
    check_equal(gv, gv_check)
    

    # TEST BOX #
    readReg = hy.ReadRegion(snapdir, 0, [*galpos-1.0, 2.0, 2.0, 2.0],
                            shape='Box', anchor='bottom',
                            verbose=True, exact=True)
    gpos = readReg.read_data('Coordinates')
    gids = readReg.read_data('ParticleIDs')
    grelpos = gpos - galpos[None, :]
    ind_box = np.nonzero(np.max(np.abs(grelpos), axis = 1) <= 1.0)[0]

    if len(ind_box) == len(grelpos):
        print("Passed")
    else:
        set_trace()

    mtot = readReg.total_in_region('Mass')

    ind_box_full = np.nonzero((np.abs(gpos_full[:, 0]-galpos[0]) <= 1.0) &
                              (np.abs(gpos_full[:, 1]-galpos[1]) <= 1.0) &
                              (np.abs(gpos_full[:, 2]-galpos[2]) <= 1.0))[0]
    print("N_rr = {:d}, N_full = {:d}"
          .format(len(gids), len(ind_sphere_full)))
    check_equal(np.sort(gids), np.sort(gids_full[ind_box_full]))


def check_equal(arr1, arr2):
    if np.all(np.equal(arr1, arr2)):
        print("Passed")
    else:
        set_trace()
        
if __name__ == "__main__":
    test_read_region()
    
