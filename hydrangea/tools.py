"""A collection of relatively simple routines for Hydrangea analysis."""

import glob
import h5py as h5
import numpy as np
import os
import hydrangea.hdf5 as hd
import hydrangea.units as hu
import ctypes as c

from astropy.io import ascii
from astropy.cosmology import Planck13
import astropy.units as u
from scipy.interpolate import interp1d

from pdb import set_trace


def get_snepshot_indices(rundir, snep_list='basic'):
    """
    Extract type, number, and aexp for snepshots from a specified list.

    Parameters
    ----------
    rundir : str
        The base directory of the simulation.
    snep_list : str, optional
        The snepshot list to extract, without extension (default: 'basic')

    Returns
    -------
    root_index : np.array(int)
        The root indices of all snepshots in specified list.
    aexp : np.array(float):
        The expansion factor of all snepshots in specified list.
    source_type : np.array(string)
        For each snepshot, whether it is a 'snap' or 'snip'.
    source_num : np.array(int)
        Indices of snepshot in its category (i.e. snap/snipshot).

    Example
    -------
    >>> import hydrangea as hy
    >>> sim = hy.objects.Simulation(index=0)
    >>> get_snepshot_indices(sim.run_dir, snep_list='z0_only')
    (array([206]), array(1.]), array(['snap'], dtype='<U4'), array([29]))
    """
    snepdir = rundir + '/sneplists/'
    file_name = snepdir + snep_list + '.dat'

    data = ascii.read(file_name)

    root_index = np.array(data['rootIndex'])
    aexp = np.array(data['aexp'])
    source_type = np.array(data['sourceType'])
    source_num = np.array(data['sourceNum'])

    return root_index, aexp, source_type, source_num


def snep_times(time_type='aexp', snep_list='allsnaps'):
    """Return the times of all snepshots in a particular list.

    By default, the expansion factors of the 30 snapshots are returned.

    Parameters
    ----------
    time_type : string
        Specified the 'time' flavour to be returned. Options are:
            'aexp' [default]: Expansion factor
            'zred': Redshift
            'age': Age of the Universe [Gyr]
            'lbt': Lookback time from z = 0 [Gyr]
    snep_list : string
        The snepshot list for which to load times. Options are:
            'z0_only': 1 snapshot at z = 0
            'regsnaps': 28 regular snapshots
            'allsnaps' [default]: 30 snapshots
            'basic': 109 snepshots at Delta_t = 125 Myr
            'default_long': basic + filler snipshots to Delta_t = 25 Myr
            'short_movie': 217 snepshots at Delta_t = 62.5 Myr
            'full_movie': 1081 snepshots at Delta_t = 12.5 Myr

    Returns
    -------
    np.array(float)
        The desired snepshot times.

    Note
    ----
    This function retrieves the target time of the snepshot, which may
    deviate slightly from the the actual output. Use read_snepshot_time()
    for the latter. It can be called directly from this library, without
    any downloaded data.

    Example
    -------
    >>> import hydrangea as hy
    >>> hy.snep_times(time_type='zred', snep_list='allsnaps')
    array([14.00321069,  6.77181334,  4.61422735,  3.51234012,  2.82506304,
        2.34762658,  1.99262408,  1.71594403,  1.49271753,  1.30777929,
        1.15130655,  1.01663116,  0.89906146,  0.79519742,  0.70250004,
        0.61903886,  0.54331181,  0.47414852,  0.410597  ,  0.36566854,
        0.35189723,  0.29742561,  0.24666517,  0.19918125,  0.15461335,
        0.11265709,  0.10063854,  0.07304807,  0.03555967,  0.        ])

    """
    curr_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
    snaptimes_file = curr_dir + 'OutputLists/' + snep_list + '.dat'
    aexp = np.array(ascii.read(snaptimes_file, format='no_header')['col1'])

    return aexp_to_time(aexp, time_type)


def aexp_to_time(aexp, time_type='age'):
    """Convert expansion factor to another time flavour.

    Parameters
    ----------
    aexp : float or np.array(float)
        (Array of) expansion factors to convert.
    time_type : string
        Specified the 'time' flavour to be returned. Options are:
            'aexp': Expansion factor
            'zred': Redshift
            'age' [default]: Age of the Universe [Gyr]
            'lbt': Lookback time from z = 0 [Gyr]

    Returns
    -------
    time : float or np.array(float)
        The desired time stamp corresponding to the input aExp(s).

    Note
    ----
    Conversions to age and lookback time assume a Planck 13 cosmology.
    For these, large input lists (> 1000 elements) are transformed via
    an interpolation from a fine-spaced grid in aexp, for speed gains.

    Example
    -------
    >>> import hydrangea as hy
    >>> hy.aexp_to_time(0.5)
    5.863165023566443
    """
    if time_type == 'aexp':
        return aexp
    if time_type == 'zred':
        return 1/aexp - 1

    if np.isscalar(aexp) or len(aexp) < 1000:
        if time_type == 'age':
            return Planck13.age(1/aexp - 1).to(u.Gyr).value
        if time_type == 'lbt':
            return Planck13.lookback_time(1/aexp - 1).to(u.Gyr).value
    else:
        aexps = np.linspace(aexp.min(), aexp.max(), 200)
        times = np.zeros(200)

        if time_type == 'age':
            times = Planck13.age(1/aexps-1).to(u.Gyr).value
        elif time_type == 'lbt':
            times = Planck13.lookback_time(1/aexps - 1).to(u.Gyr).value
        else:
            print("I do not understand time_type '{:s}'!"
                  .format(time_type))
            set_trace()

        csi_time = interp1d(aexps, times, kind='cubic')
        return csi_time(aexp)


def get_m_dm(file_name, units='astro'):
    """Retrieve the DM particle mass from a particle file.

    Parameters
    ----------
    file_name : str
        The full name of the file from which to look up the mass.
    units : str, optional
        The desired unit system for the result (default: 'astro', i.e. M_sun)

    Returns
    -------
    m_dm : float
        The mass of the DM particle.

    Note
    ----
    This value is retrievable from any snapshot file,
    and is identical across different outputs from the same simulation.
    It does, however, vary slightly between different simulations.

    Example
    -------
    >>> import hydrangea as hy
    >>> sim = hy.objects.Simulation(index=0)
    >>> hy.get_m_dm(sim.get_snap_file(29))
    0.0009615178624612349
    """
    m_dm = hd.read_attribute(file_name, 'Header', 'MassTable')[1]
    if units == 'data':
        return m_dm
    clean_factor = 1.0/hd.read_attribute(file_name, 'Header', 'HubbleParam')
    m_dm *= clean_factor
    if units in ['clean', 'astro']:
        return m_dm
    m_dm *= 1e10 * hu.SOLAR_MASS
    if units == 'si':
        return m_dm
    if units == 'cgs':
        return m_dm * 1e3
    else:
        print(f"Unsupported unit system '{units}'!")


def get_m_baryon(file_name, units='astro'):
    """Retrieve the initial baryon mass from a particle file."""
    m_dm = get_m_dm(file_name, units=units)
    omega_matter = hd.read_attribute(file_name, 'Header', 'Omega0')
    omega_baryon = hd.read_attribute(file_name, 'Header', 'Omega_Baryon')
    return m_dm * omega_baryon/omega_matter / (1-omega_baryon/omega_matter)


def ind_to_block(indices, offsets, lengths=None):
    """Find the block for a set of indices in an offset-separated list.

    Each block i contains elements from indices offsets[i] up to and
    including offsets[i]+length[i]-1. This can, for example, be used to
    translate particle indices in a Subfind-ID list into the corresponding
    subhalo or FOF index.

    Parameters
    ----------
    indices : ndarray(int)
        The indices of elements whose blocks to find.
    offsets : ndarray(int)
        The index of the first element in each block.
    lengths : ndarray(int), optional
        The number of elements in each block. If None (default), it is
        assumed that the blocks are contiguous, i.e. all elements between
        offsets[i] and offsets[i+1] belong to block i.

    Returns
    -------
    blocks : ndarray(int)
        The block index for each input element (-1 if not found).

    Note
    ----
    When lengths is not provided, it is advisable to append a `coda'
    to offsets, i.e. a trailing entry with the total number of elements
    assigned to blocks. This enables the correct identification of
    (potential) input elements beyond the range of the last block.

    Example
    -------
    >>> import hydrangea as hy
    >>> indices = np.array([1, 5, 6, 79, 81])
    >>> offsets = np.array([0, 6, 10, 50, 80, 100])
    >>> lengths = np.array([6, 4, 20, 20, 2])
    >>> hy.ind_to_block(indices, offsets, lengths)
    array([ 0, 0,  1, -1,  4], dtype=int32)
    """
    block_guess = np.searchsorted(offsets, indices, side='right')-1

    if lengths is None:
        ind_good_guess = np.nonzero(block_guess >= 0)[0]
    else:
        ind_good_guess = np.nonzero(
            (indices >= offsets[block_guess]) &
            (indices < offsets[block_guess] + lengths[block_guess]))[0]

    block_index = np.zeros(len(indices), dtype=np.int32) - 1
    block_index[ind_good_guess] = block_guess[ind_good_guess]

    return block_index


def form_files(sim_dir, index, types='sub', snep_type='snap'):
    """Create the file names for different output types and snapshots.

    This is a convenience function to avoid having to manually construct
    the HDF5 file names. Works for both Hydrangea and other Eagle simulations.

    Parameters
    ----------
    sim_dir : string
        The simulation base directory.
    index : int
        The snap- or snipshot index for which to create file names.
    types : string or list of strings, optional
        The type(s) of output for which to construct file names. Can be
        one or more of the following (separate by space if multiple):
            'sub'     --> Subfind subhalo catalogue (default)
            'snap'    --> Raw snapshot
            'fof'     --> Subfind FOF catalogue
            'subpart' --> Subfind particle catalogue
            'partmag' --> Magnitudes of star particles (Hydrangea only)
    snep_type : string, optional
        Snepshot type to construct file(s) for: 'snap[shot]' (default)
        or 'snip[shot]' (both short and long forms accepted, capitalization
        ignored).

    Returns
    -------
    file(s) : string or None (or list thereof)
        Name(s) of the *first* HDF5 file of the desired output(s). If
        multiple types are specified, the filenames are in the same order
        as the input. None is returned for file names that could not be found.
    """
    # Consistency checks on input parameters:
    s_type = snep_type.lower()[:4]
    if s_type not in ['snip', 'snap']:
        print("Invalid snep_type '{:s}'!" .format(snep_type))
        set_trace()

    # Split possibly multiple type inputs into list
    types_list = types.split()
    files_list = []
    z_string = None

    # Snapshot index as a triple-zero-padded string:
    index_string = '{:03d}' .format(index)

    # Translate each type into a directory and file prefix:
    for itype in types_list:
        if itype.lower() == 'snap':
            if s_type == 'snap':
                names = ('snapshot', 'snap')
            else:
                names = ('snipshot', 'snip')

        elif itype.lower() == 'sub':
            names = ('groups', 'eagle_subfind_tab')
        elif itype.lower() == 'subpart':
            names = ('particledata', 'eagle_subfind_particles')
        elif itype.lower() == 'fof':
            names = ('groups', 'group_tab')
        elif itype.lower() == 'partmag':
            names = ('snapshot', 'partMags_EMILES_PDXX_DUST_CH')
            if s_type != 'snap':
                raise Exception("Stellar magnitudes only available for "
                                "snapshots.")
        else:
            raise Exception("Data type '" + itype + "' is not understood."
                            "Please try another one.")

        if z_string is None:
            z_string = _find_z_string(sim_dir, names[0], index_string)
            if z_string is None:
                print("Error determining z-string...")
                set_trace()

        # Build appropriate file name (different for magnitudes):
        if itype == 'partmag':
            ifile = ("{0:s}/data/stars_extra/{1:s}_{2:s}_{3:s}/"
                     "{4:s}_{2:s}_{3:s}.0.hdf5"
                     .format(sim_dir,
                             names[0], index_string, z_string, names[1]))
        else:
            ifile = ("{0:s}/data/{1:s}_{2:s}_{3:s}/{4:s}_{2:s}_{3:s}.0.hdf5"
                     .format(sim_dir,
                             names[0], index_string, z_string, names[1]))

        files_list.append(ifile)

    # If input is single string, output should be too
    if len(types_list) == 1:
        files_list = files_list[0]

    return files_list


def _find_z_string(sim_dir, dir_type, index_string):
    """Determine the redshift string for a particular output."""
    s_dir = glob.glob(sim_dir + 'data/' + dir_type + '_' + index_string + '_*')
    if s_dir:
        dir_name = (s_dir[0].split('/'))[-1]
        z_string = (dir_name.split('_'))[-1]
        return z_string
    else:
        return None


def sum_bins(*args, **kwargs):
    """Sum up a quantity split by multiple indices.

    Parameters
    ----------
    quant : ndarray
        The quantity to sum up (float32)
    indices : ndarrays (1-8)
        The index or indices according to which to sum the quantity.

    Returns
    -------
    sum_array : ndarray (float32)
        The sum of the quantity in each combination of supplied indices.

    Note
    ----
    The numpy.histogramdd provides similar functionality, but this function
    is typically somewhat faster. It uses the Kahan summation algorithm
    (https://en.wikipedia.org/wiki/Kahan_summation_algorithm) for
    accurate float32 summation even when summing many elements.
    """
    quant = args[0]
    num_indices = len(args) - 1
    if num_indices < 1:
        print("You need to provide at least one index!")
        set_trace()
    if num_indices > 8:
        print("Can currently only bin by up to 8 indices, not {num_indices}.")
        set_trace()

    num_elem = len(quant)

    if 'max_indices' in kwargs:
        max_indices = kwargs['max_indices']
        if len(max_indices) != num_indices:
            print("Length of max_indices does not match indices...")
            set_trace()
    else:
        max_indices = np.zeros(num_indices, dtype=int)
        for ii in range(num_indices):
            max_indices[ii] = np.max(args[ii+1])

    quant_binned = np.zeros((max_indices + 1), dtype=np.float32)
    kahan_temp = np.zeros((max_indices + 1), dtype=np.float32)

    # Convert metadata to C-compatible format
    c_num_elem = c.c_long(num_elem)
    c_nbins_all = [c.c_int(max_indices[ii]+1) for ii in range(num_indices)]
    for ii in range(num_indices, 8):
        c_nbins_all.append(c.c_int(1))

    # Get pointers to in- and output arrays
    c_nbins_pointers = [c.addressof(_c_nbins) for _c_nbins in c_nbins_all]

    if quant.dtype != np.float32:
        print(id(quant))
        quant = np.array(quant, dtype=np.float32)
    quant_p = quant.ctypes.data_as(c.c_void_p)

    indices_p = []
    transargs = [None]*8
    dummy_array = np.zeros(num_elem, dtype=np.int32)
    for ii in range(8):
        if ii < num_indices:
            if args[ii+1].dtype == np.int32:
                indices_p.append(args[ii+1].ctypes.data_as(c.c_void_p))
            else:
                transargs[ii] = args[ii+1].astype(np.int32)
                indices_p.append(transargs[ii].ctypes.data_as(c.c_void_p))
        else:
            indices_p.append(dummy_array.ctypes.data_as(c.c_void_p))
    result_p = quant_binned.ctypes.data_as(c.c_void_p)
    kahan_p = kahan_temp.ctypes.data_as(c.c_void_p)

    nargs = 20
    myargv = c.c_void_p * nargs
    argv = myargv(c.addressof(c_num_elem),
                  *c_nbins_pointers,
                  quant_p,
                  *indices_p,
                  result_p, kahan_p)

    object_dir = os.path.dirname(os.path.realpath(__file__)) + "/clib/"
    lib = c.cdll.LoadLibrary(object_dir + 'sumbins.so')
    lib.sumbins(nargs, argv)
    return quant_binned
