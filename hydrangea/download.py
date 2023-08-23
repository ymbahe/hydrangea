# -*- coding: utf-8 -*-
"""Module to download files from the MPCDF virgo archive."""

import requests
import os
import time
from pdb import set_trace
import numpy as np
import hydrangea as hy

IGNORE_EXISTING = False
BASE_DIR = None


def download_data(base_dir, suite='10r200', sim=None, index=None,
                  types=None, sim_type='HYDRO', reload=False):
    """Download a specified collection of data.

    This is the main public function to retrieve the data. Users can
    specify the range of simulations, output indices, and data
    products to download.

    Warning
    -------
    With this function, it is possible to download huge amounts of data
    (up to several hundred TB). Users should ensure that there is
    sufficient disk space and bandwidth to accommodate this download,
    and that unneccessary data are not downloaded by accident.

    Parameters
    ----------
    base_dir : str
        Local base directory of the data repository. The data will
        be downloaded to this location, in a structure mirrored from
        the central MPCDF archive.
    suite : str, optional
        Simulation suite to download data from. Options are ``'10r200'``
        (default), ``'5r200'``, ``'10r200_only'``, and ``'5r200_only'``.
        The ..._only options do not default to loading simulation data
        from the other suite if the specified simulation is not
        included in the specified one.
    sim : int or list/array of ints or ``None``, optional
        The index (or indices) of the simulation(s) to download. If
        ``None`` (default), no data is downloaded (this is intended
        as a precaution against accidentally downloading an unwanted
        simulation).
    index : int or list/array of ints or ``None``, optional
        For output-specific data (particle or subfind catalogues),
        specify the output index (or indices) to download. If ``None``
        (default), no data is downloaded in this case (intended as a
        precaution against accidentally downloading unwanted simulation
        outputs).
    types : str or ``None``, optional
        Types of data to download. Multiple types can be specified with
        an empty space between them. Options are:

        - ``'all'`` : All data (equivalent to ``'highlev logs magnitudes
          snapshot subfind'``)
        - ``'evolution_tables'`` : Pre-compiled tables of galaxy properties
          in all snapshots
        - ``'highlev'`` : All high-level simulation galaxy catalogues
          (equivalent to ``'spiderweb subhalo_extra evolution_tables
          galaxy_coordinates galaxy_paths'``)
        - ``'galaxy_coordinates'`` : Pre-compiled tables of galaxy positions
          and velocities
        - ``'galaxy_paths'`` : Approximate galaxy positions and velocities
          in all outputs, including snipshots
        - ``'logs'`` : Log, code, and parameter files
        - ``'magnitudes'`` : Stellar particle magnitude catalogue(s)
        - ``'subhalo_magnitudes'``: Subhalo integrated magnitude catalogue(s)
        - ``'snapshot'`` : Snapshot particle catalogue(s)
        - ``'snipshot'`` : Snipshot particle catalogue(s)
          (*not yet implemented*)
        - ``'spiderweb'`` : Spiderweb tables to trace galaxies between
          snapshots
        - ``'subfind'`` : Subfind catalogue(s)
        - ``'subhalo_extra'`` : Additional subhalo properties that are
          not stored in the subfind catalogue
    sim_type : str, optional
        Type of simulation to download data from. Options are
        'HYDRO' (default, for the simulations with baryons) or 'DM'
        for the corresponding gravity-only simulations (case-insensitive).
    reload : bool, optional
        Re-download already existing files (default: ``False``).
    """
    global IGNORE_EXISTING
    global BASE_DIR

    if not reload:
        IGNORE_EXISTING = True

    sim_type = sim_type.upper()

    if not base_dir.endswith('/'):
        base_dir = base_dir + '/'
    BASE_DIR = base_dir

    # Make sure request is sensible and complete
    if sim is None or types is None:
        print("Please specify what to download:")
        print("sim   --> Simulation (range)")
        print("types --> Data types")
        return

    # Resolve shorthands for collections:
    types = types.replace('all', ('logs snapshot subfind magnitudes '
                                  'subhalo_magnitudes highlev'))
    types = types.replace('highlev', ('spiderweb subhalo_extra match '
                                      'evolution_tables galaxy_coordinates '
                                      'galaxy_paths'))

    type_list = [_t.lower() for _t in types.split()]
    if index is None:
        for itype in type_list:
            if itype in ['snapshot', 'subfind', 'magnitudes',
                         'subhalo_magnitudes', 'cantor']:
                print(f"To download data of type '{itype}', you must "
                      "specify which output indices to retrieve...")
                return

    sim = _to_array(sim)
    index = _to_array(index)

    if index is not None:
        zstrings = get_zstrings(index)
    if np.isscalar(index):
        zstrings = [zstrings]

    for isim in sim:
        sim_base = f'{suite}/CE-{isim}/{sim_type}'
        for itype in type_list:
            if itype == 'snapshot':
                for iindex, index_curr in enumerate(index):
                    download_snapshot(sim_base, index_curr, zstrings[iindex])
            elif itype == 'snipshot':
                for iindex, index_curr in enumerate(index):
                    download_snipshot(sim_base, index_curr, zstrings[iindex])
            elif itype == 'subfind':
                for iindex, index_curr in enumerate(index):
                    download_subfind(sim_base, index_curr, zstrings[iindex])
            elif itype == 'magnitudes':
                for iindex, index_curr in enumerate(index):
                    download_magnitudes(sim_base, index_curr,
                                        zstrings[iindex])
            elif itype == 'subhalo_magnitudes':
                for iindex, index_curr in enumerate(index):
                    download_subhalo_magnitudes(sim_base, index_curr,
                                                zstrings[iindex])
            elif itype == 'spiderweb':
                download_spiderweb(sim_base, target='tables')
            elif itype == 'spiderweb_links':
                download_spiderweb(sim_base, target='links')
            elif itype == 'run_logs':
                download_run_logs(sim_base)
            elif itype == 'step_logs':
                download_step_logs(sim_base)

            elif itype == 'subhalo_extra':
                download_highlev(sim_base, 'SubhaloExtra.hdf5')                   
            elif itype == 'evolution_tables':
                download_highlev(sim_base, 'FullGalaxyTables.hdf5')
            elif itype == 'galaxy_coordinates':
                download_highlev(sim_base, 'GalaxyPositionsSnap.hdf5')
            elif itype == 'galaxy_paths':
                download_highlev(sim_base, 'GalaxyPaths.hdf5')


def download_snapshot(sim_base, index, zstring):
    """Wrapper to download a snapshot file set."""
    file_name = (f'{sim_base}/data/snapshot_{index:03d}_z{zstring}/'
                 f'snap_{index:03d}_z{zstring}.0.hdf5')
    map_name = (f'{sim_base}/data/snapshot_{index:03d}_z{zstring}/'
                'ParticleMap.hdf5')

    download_files(file_name)
    download_single_file(map_name)


def download_snipshot(sim_base, index, zstring):
    """Wrapper to download a snipshot file set."""
    file_name = (f'{sim_base}/data/snipshot_{index:03d}_z{zstring}/'
                 f'snip_{index:03d}_z{zstring}.0.hdf5')
    map_name = (f'{sim_base}/data/snipshot_{index:03d}_z{zstring}/'
                'ParticleMap.hdf5')

    download_files(file_name)
    download_single_file(map_name)

def download_subfind(sim_base, index, zstring):
    """Wrapper to download a subfind file set."""
    print(f"--- Downloading subfind catalogue {index} (sim {sim_base}) ---")
    file_name = (f'{sim_base}/data/groups_{index:03d}_z{zstring}/'
                 f'eagle_subfind_tab_{index:03d}_z{zstring}.0.hdf5')
    download_files(file_name)


def download_magnitudes(sim_base, index, zstring):
    """Wrapper to download a magnitude file set."""
    print(f"--- Downloading stellar magnitudes {index} (sim {sim_base}) ---")
    file_name = (f'{sim_base}/data/stars_extra/'
                 f'snapshot_{index:03d}_z{zstring}/'
                 f'partMags_EMILES_PDXX_DUST_CH_{index:03d}_z{zstring}.0.hdf5')
    download_files(file_name)

def download_subhalo_magnitudes(sim_base, index, zstring):
    """Wrapper to download a subhalo magnitude file set."""
    print(f"--- Downloading subhalo magnitudes {index} (sim {sim_base}) ---")
    file_name = (f'{sim_base}/data/stars_extra/'
                 f'groups_{index:03d}_z{zstring}/'
                 f'galaxyMagnitudes_EMILES_PDXX_DUST_CH_{index:03d}_'
                 f'z{zstring}.0.hdf5')
    download_files(file_name)

def download_spiderweb(sim_base, target='tables'):
    """Download spiderweb data."""
    print(f"--- Downloading SpiderWeb catalogue (sim {sim_base}) ---")
    if target == 'tables':
        cat_file = f'{sim_base}/highlev/SpiderwebTables.hdf5'
    elif target == 'links':
        cat_file = f'{sim_base}/highlev/SpiderwebLinks.hdf5'
    else:
        print(f"Spiderweb sub-type '{target}' not understood...")
        return

    download_single_file(cat_file)


def download_run_logs(sim_base):
    """Download log files."""
    print(f"--- Downloading log files (sim {sim_base}) ---")
    download_single_file(f'{sim_base}/logfiles.tar.gz')
    download_single_file(f'{sim_base}/logfiles.toc')


def download_step_logs(sim_base):
    """Download time step logs."""
    print(f"--- Downloading time step log files (sim {sim_base}) ---")

    txt_list = ['balance.txt', 'blackholes.txt', 'cpu.txt',
                'eagle_blackholes.txt', 'energy.txt', 'info.txt',
                'metals_gas.txt', 'metals_sf.txt', 'metals_stars.txt',
                'metals_tot.txt', 'sfr.txt', 'SNIa.txt',
                'StellarEvolutionAGB.txt', 'StellarEvolutionIa.txt',
                'StellarEvolutionII.txt', 'StellarEvolutionTotal.txt',
                'stellar_feedback.txt', 'timebins.txt', 'timings.txt']
    for itxt in txt_list:
        download_single_file(f'{sim_base}/data/{itxt}')


def download_highlev(sim_base, target="FullGalaxyTables.hdf5"):
    """Download a specified file from highlev directory."""
    print(f" --- Downloading highlev/{target} (sim {sim_base}) ---")
    download_single_file(f'{sim_base}/highlev/{target}')


def download_files(file_name):
    """Download a set of related files."""
    download_single_file(file_name)
    nfiles = hy.hdf5.read_attribute(BASE_DIR + file_name,
                                    'Header', 'NumFilesPerSnapshot')
    if nfiles is None:
        print("Error: could not determine number of files in set...")
        set_trace()
    for ifile in range(1, nfiles):
        curr_file_name = _swap_file_name(file_name, ifile)
        download_single_file(curr_file_name)


def get_zstrings(indices, sneplist='allsnaps'):
    """Generate the redshift strings for specified output(s)."""
    zred = hy.snep_times('zred', sneplist)[indices]
    if np.isscalar(zred):
        return f'{zred:07.3f}'.replace('.', 'p')
    else:
        return [f'{_z:07.3f}'.replace('.', 'p') for _z in zred]


def download_single_file(source):
    """Download one single file to a specified location.

    This is the low-level 'worker' function that is called by the
    wrapper 'download_files()'.
    """
    # Save file to appropriate place

    save_loc = BASE_DIR + source
    if not os.path.isdir(os.path.dirname(save_loc)):
        os.makedirs(os.path.dirname(save_loc))

    if IGNORE_EXISTING:
        if os.path.exists(save_loc):
            print(f"File {source} already exists -- skipping.")
            return

    stime = time.time()

    # Get file from archive
    print(f"Downloading file {source}...", end='')
    url = f'https://hydrangea.mpcdf.mpg.de/{source}'
    resp = requests.head(url, auth=('public', 'flower'))
    size_bytes = int(resp.headers['Content-Length'])
    size_mb = size_bytes / (1024*1024)
    print(f' [{size_mb:.3f} MB]')
    r = requests.get(url, auth=('public', 'flower'), allow_redirects=True,
                     stream=True)
    r.raise_for_status()
    if r.status_code > 400:
        print("Looks like something went wrong with the download of "
              f"'{source}'...")

    with open(save_loc, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:   # filter out keep-alive new chunks
                f.write(chunk)

    duration = time.time() - stime
    rate = size_mb/duration
    print(f'   ---> {duration:1f} sec., {rate:.3f} MB/s')


def _to_array(x):
    if np.isscalar(x):
        return np.array([x])
    if isinstance(x, tuple):
        return np.arange(x[0], x[1]+1)
    return x


def _swap_file_name(base_name, number):
    """Form sequentially numbered file names (convenience function)."""
    name_parts = base_name.split('.')
    name_parts[-2] = str(number)
    new_name = '.'.join(name_parts)
    return new_name


if __name__ == '__main__':
    """TEMPORARY section to test this module."""
    download_single_file('10r200/CE-0/HYDRO/data/snapshot_029_z000p000/',
                         '/Users/bahe/publictest/')
