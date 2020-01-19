#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module to download files from the virgo archive."""

import requests
import os
import time
from pdb import set_trace
import numpy as np
import hydrangea as hy

IGNORE_EXISTING = False
BASE_DIR = None

def download_data(base_dir, suite='10r200', sim=None, index=None,
                  types=None, index_type='allsnaps', sim_type='HYDRO',
                  reload=False):
    """Download a specified collection of data.
    
    This is the main public function to retrieve the data. Users can
    specify the range of simulations, output indices, and data 
    products to download.
    """

    global IGNORE_EXISTING
    global BASE_DIR

    if not reload:
        IGNORE_EXISTING = True

    if not base_dir.endswith('/'): base_dir = base_dir + '/'
    BASE_DIR = base_dir

    # Make sure request is sensible and complete
    if sim is None or types is None:
        print("Please specify what to download:")
        print("sim   --> Simulation (range)")
        print("types --> Data types")
        return    

    # Resolve shorthands for collections:
    types = types.replace('all', ('logs snapshot subfind magnitudes highlev'))
    types = types.replace('highlev', ('spiderweb subhalo_extra match '
                                      'evolution_tables galaxy_coordinates '
                                      'galaxy_paths'))

    type_list = [_t.lower() for _t in types.split()]
    if index is None:
        for itype in type_list:
            if itype in ['snapshot', 'subfind', 'magnitudes',
                         'cantor']:
                print(f"To download data of type '{itype}', you must "
                      "specify which output indices to retrieve...")
    
    sim = _to_array(sim)
    index = _to_array(index)

    if index is not None:
        zstrings = get_zstrings(index)
    
    for isim in sim:
        sim_base = f'{suite}/CE-{isim}/{sim_type}'
        for itype in type_list:
            if itype == 'snapshot':
                for iindex in index:
                    download_snapshot(base_dir, sim_base, 
                                      iindex, zstrings[iindex])
            elif itype == 'subfind':
                for iindex in index:
                    download_subfind(base_dir, sim_base,
                                     iindex, zstrings[iindex])
            elif itype == 'magnitudes':
                for iindex in index:
                    download_magnitudes(base_dir, sim_base,
                                        iindex, zstrings[iindex])
            elif itype == 'spiderweb':
                download_spiderweb(base_dir, sim_base,
                                   target='tables')
            elif itype == 'spiderweb_links':
                download_spiderweb(base_dir, sim_base,
                                   target='links')
            elif itype == 'logs':
                download_logs(base_dir, sim_base)

            elif itype == 'subhalo_extra':
                download_highlev(base_dir, sim_base, 'SubhaloExtra.hdf5')
            elif itype == 'match':
                if sim_type == 'HYDRO':
                    download_highlev(base_dir, sim_base, 'MatchInDM.hdf5')
                elif sim_type == 'DM':
                    download_highlev(base_dir, sim_base, 'MatchInHydro.hdf5')
                else:
                    print(f"Unknown simulation type {sim_type}!")
                    set_trace()
            elif itype == 'evolution_tables':
                download_highlev(base_dir, sim_base, 'FullGalaxyTables.hdf5')
            elif itype == 'galaxy_coordinates':
                download_highlev(base_dir, sim_base, 
                                 'GalaxyPositionsSnap.hdf5')
            elif itype == 'galaxy_paths':
                download_highlev(base_dir, sim_base, 
                                 'GalaxyPaths.hdf5')
                

def download_snapshot(base_dir, sim_base, index, zstring):
    """Wrapper to download a snapshot file set."""
    file_name = (f'{sim_base}/data/snapshot_{index:03d}_z{zstring}/'
                 f'snap_{index:03d}_z{zstring}.0.hdf5')
    map_name = (f'{sim_base}/data/snapshot_{index:03d}_z{zstring}/'
                'ParticleMap.hdf5')

    download_files(file_name, base_dir)
    download_single_file(map_name, base_dir)


def download_subfind(base_dir, sim_base, index, zstring):
    """Wrapper to download a subfind file set."""
    print(f"--- Downloading subfind catalogue {index} (sim {sim_base}) ---")
    file_name = (f'{sim_base}/data/groups_{index:03d}_z{zstring}/'
                 f'eagle_subfind_tab_{index:03d}_z{zstring}.0.hdf5')
    download_files(file_name, base_dir)


def download_magnitudes(base_dir, sim_base, index, zstring):
    """Wrapper to download a magnitude file set."""
    print(f"--- Downloading stellar magnitudes {index} (sim {sim_base}) ---")
    file_name = (f'{sim_base}/data/stars_extra/'
                 f'snapshot_{index:03d}_z{zstring}/'
                 f'partMags_EMILES_PDXX_DUST_CH_{index:03d}_z{zstring}.0.hdf5')
    download_files(file_name, base_dir)


def download_spiderweb(base_dir, sim_base, target='tables'):
    """Download spiderweb data"""
    print(f"--- Downloading SpiderWeb catalogue (sim {sim_base}) ---")
    if target == 'tables':
        cat_file = f'{sim_base}/highlev/SpiderwebTables.hdf5'
    elif target == 'links':
        cat_file = f'{sim_base}/highlev/SpiderwebLinks.hdf5'
    else:
        print(f"Spiderweb sub-type '{target}' not understood...")
        return

    download_single_file(cat_file, base_dir)
    
def download_logs(base_dir, sim_base):
    """Download log files."""
    print(f"--- Downloading log files (sim {sim_base}) ---")
    download_single_file(f'{sim_base}/logfiles.tar.gz', base_dir)
    download_single_file(f'{sim_base}/logfiles.toc', base_dir)
    download_single_file(f'{sim_base}/logfiles.toc', base_dir)

    txt_list = ['balance.txt', 'blackholes.txt', 'cpu.txt', 
                'eagle_blackholes.txt', 'energy.txt', 'info.txt',
                'metals_gas.txt', 'metals_sf.txt', 'metals_stars.txt',
                'metals_tot.txt', 'sfr.txt', 'SNIa.txt',
                'StellarEvolutionAGB.txt', 'StellarEvolutionIa.txt',
                'StellarEvolutionII.txt', 'StellarEvolutionTotal.txt',
                'stellar_feedback.txt', 'timebins.txt', 'timings.txt']
    for itxt in txt_list:
        download_single_file(f'{sim_base}/data/{itxt}', base_dir)


def download_highlev(base_dir, sim_base, target="FullGalaxyTables.hdf5"):
    """Download a specified file from highlev directory."""
    print(f" --- Downloading highlev/{target} (sim {sim_base}) ---")
    download_single_file(f'{sim_base}/highlev/{target}')


def download_files(file_name, base_dir):
    """Download a set of related files."""
    download_single_file(file_name, base_dir)
    nfiles = hy.hdf5.read_attribute(base_dir + file_name,
                                    'Header', 'NumFilesPerSnapshot')
    if nfiles is None:
        print("Error: could not determine number of files in set...")
        set_trace()
    for ifile in range(1, nfiles):
        curr_file_name = _swap_file_name(file_name, ifile)
        download_single_file(curr_file_name, base_dir)


def get_zstrings(indices, sneplist='allsnaps'):
    """Generate the redshift strings for specified output(s)."""
    zred = hy.snep_times('zred', sneplist)[indices]
    if np.isscalar(zred):
        return f'{zred:07.3f}'.replace('.', 'p')
    else:
        return [f'{_z:07.3f}'.replace('.', 'p') for _z in zred]

def download_single_file(source, base_dir):
    """Download one single file to a specified location.
    
    This is the low-level 'worker' function that is called by the
    wrapper 'download_files()'.
    """ 
    
    if IGNORE_EXISTING:
        if os.path.exists(BASE_DIR + source):
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
    r = requests.get(url, auth=('public', 'flower'), allow_redirects=True)
    if r.status_code > 400:
        print("Looks like something went wrong with the download of " 
              f"'{source}'...")

    # Save it to appropriate place
    if not base_dir.endswith('/'):
        base_dir = base_dir + '/'
    save_loc = base_dir + source
    if not os.path.isdir(os.path.dirname(save_loc)):
        os.makedirs(os.path.dirname(save_loc))
    open(BASE_DIR + source, 'wb').write(r.content)                     

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
    download_single_file('10r200/CE-0/HYDRO/data/snapshot_029_z000p000/',
                         '/Users/bahe/publictest/')
       

