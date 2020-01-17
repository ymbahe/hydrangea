#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module to download files from the virgo archive."""

import requests
import os
from pdb import set_trace
import numpy as np
import hydrangea as hy

def download_data(base_dir, suite='10r200', sim=None, index=None,
                  types=None, index_type='allsnaps', sim_type='HYDRO'):
    """Download a specified collection of data.
    
    This is the main public function to retrieve the data. Users can
    specify the range of simulations, output indices, and data 
    products to download.
    """

    # Make sure request is sensible and complete
    if sim is None or types is None:
        print("Please specify what to download:")
        print("sim   --> Simulation (range)")
        print("types --> Data types")
        return    
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
        sim_base = f'{suite}/CE-{isim}/{sim_type}/'
        for itype in type_list:
            if itype == 'snapshot':
                for iindex in index:
                    download_snapshot(base_dir, sim_base, 
                                      iindex, zstrings[iindex])
            elif itype == 'subfind':
                for iindex in index:
                    download_subfind(base_dir, sim_base,
                                     iindex, zstrings[iindex])
            elif iitype == 'magnitudes':
                for iindex in index:
                    download_magnitudes(base_dir, sim_base,
                                        iindex, zstrings[iindex])
            elif iitype == 'spiderweb':
                download_spiderweb(base_dir, sim_base,
                                   target='tables')
            elif iitype == 'spiderweb_links':
                download_spiderweb(base_dir, sim_base,
                                   target='links')
            elif iitype == 'logs':
                download_logs(base_dir, sim_base)

def download_snapshot(base_dir, sim_base, index, zstring):
    """Wrapper to download a snapshot file set."""
    file_name = (f'{sim_base}/data/snapshot_{index:03d}_z{zstring}/'
                 f'snap_{index:03d}_z{zstring}.0.hdf5')
    map_name = (f'{sim_base}/data/snapshot_{index:03d}_z{zstring}/'
                'ParticleMap.hdf5')
    download_files(file_name, base_dir, multi=True)
    download_single_file(map_name, base_dir)


def download_files(base_dir, file_name, multi=False):
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
    # Get file from archive
    r = requests.get(f'https://hydrangea.mpcdf.mpg.de/{source}',
                     auth=('public', 'flower'), allow_redirects=True)
    if r.status_code > 400:
        print("Looks like something went wrong with the download of " 
              f"'{source}'...")

    # Save it to appropriate place
    if not base_dir.endswith('/'):
        base_dir = base_dir + '/'
    save_loc = base_dir + source
    if not os.path.isdir(os.path.dirname(save_loc)):
        os.makedirs(os.path.dirname(save_loc))
    open(base_dir + '/' + source, 'wb').write(r.content)                     


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
       

