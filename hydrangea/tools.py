"""A collection of relatively simple routines for Hydrangea analysis."""

import glob
import h5py as h5
import numpy as np
import os
import hydrangea.hdf5 as hd

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
    for the latter.
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


def get_astro_conv(file_name, dataset_name):
    """Get the conversion factor to astronomical units for a data set.

    Parameters
    ----------
    file_name : str
        The file containing the target data set.
    dataset_name : str
        The data set for which to compute the conversion factor.

    Returns
    -------
    conv : float or None
        Conversion factor (None if it could not be determined).
    """
    f = h5.File(file_name, 'r')
    dSet = f[dataset_name]
    try:
        hscale_exponent = dSet.attrs["h-scale-exponent"]
        ascale_exponent = dSet.attrs["aexp-scale-exponent"]

        header = f["/Header"]
        aexp = header.attrs["ExpansionFactor"]
        h_hubble = header.attrs["HubbleParam"]
    except KeyError:
        f.close()
        return None
    f.close()

    astro_conv = aexp**ascale_exponent * h_hubble**hscale_exponent
    return astro_conv


def get_m_dm(file_name, astro=True):
    """Retrieve the DM particle mass from a particle file."""
    m_dm = hd.read_attribute(file_name, 'Header', 'MassTable')[0]
    if astro:
        m_dm /= hd.read_attribute(file_name, 'Header', 'HubbleParam')
    return m_dm


def get_m_baryon(file_name, astro=True):
    """Retrieve the initial baryon mass from a particle file."""
    m_dm = get_m_dm(file_name, astro=astro)
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
