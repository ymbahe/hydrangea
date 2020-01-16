#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Classes for efficient working with simulations and snapshots."""

from pdb import set_trace

import hydrangea as hy

BASE_DIR = '/virgo/simulations/Hydrangea/'


class Simulation:
    """Representation of one simulation."""

    def __init__(self, run_dir=None, index=None, base_dir=BASE_DIR,
                 sim_type='Hydro', suite_name='10r200'):
        """Object constructor."""
        if run_dir is None:
            if base_dir is None or index is None:
                print("You need to specify either run_dir or "
                      "all the building blocks...")
                set_trace()
            if base_dir.endswith('/'):
                base_dir = base_dir[:-1]
            self.run_dir = ('{:s}/{:s}/CE-{:d}/{:s}/'
                            .format(base_dir, suite_name, index,
                                    sim_type.upper()))
        else:
            if not run_dir.endswith('/'):
                run_dir = run_dir + '/'
            self.run_dir = run_dir

        self.high_level_dir = self.run_dir + 'highlev/'
        self.cantor_dir = self.high_level_dir + 'Cantor/'

        self.fgt_loc = self.high_level_dir + 'FullGalaxyTables.hdf5'
        self.gps_loc = self.high_level_dir + 'GalaxyPositionsSnap.hdf5'
        self.galaxy_path_loc = self.high_level_dir + 'GalaxyPaths.hdf5'
        self.spider_loc = self.high_level_dir + 'SpiderwebTables.hdf5'
        self.sh_extra_loc = self.high_level_dir + 'SubhaloExtra.hdf5'

    def get_snap_file(self, index):
        """Form the (first) file of a given snapshot."""
        return hy.form_files(self.run_dir, index, types='snap')

    def get_subfind_file(self, index):
        """Form the (first) file of a given subfind catalogue."""
        return hy.form_files(self.run_dir, index, types='sub')