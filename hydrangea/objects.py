#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Classes for efficient working with simulations and snapshots."""

from pdb import set_trace
import os
import numpy as np

import hydrangea as hy

BASE_DIR = hy.local.BASE_DIR


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

    def get_snapshot_file(self, index):
        """Form the (first) file of a given snapshot."""
        return hy.form_files(self.run_dir, index, types='snap')

    def get_subfind_file(self, index):
        """Form the (first) file of a given subfind catalogue."""
        return hy.form_files(self.run_dir, index, types='sub')

    def get_snipshot_file(self, index):
        """Form the (first) file of a given snipshot."""
        return hy.form_files(self.run_dir, index, types='snap',
                             snep_type='snip')

    def get_snap_file(self, index):
        """Shorthand alias for get_snapshot_file()."""
        self.get_snapshot_file(self, index)

    def get_snip_file(self, index):
        """Shorthand alias for get_snipshot_file()."""
        self.get_snipshot_file(self, index)


class Cantor:
    """Representation of a Cantor catalogue.

    This object can be used to conveniently load and store the data
    from a Cantor catalogue.

    Similar to SplitFile, it provides the
    option to specify a reading range or explicit indices of which
    objects to read data for. Both provide the same behaviour for the
    'main' catalogue. However, for the 'extra' catalogue, indices
    are translated to the appropriate extra index, whereas ranges
    (and full catalogues) directly load the entries in the extra
    catalogue as well.
    """

    def __init__(self, sim, index=None, snap_type='snap', cantor_dir=None,
                 group='Subhalo', read_range=None, read_index=None,
                 index_dim=0):
        """Object constructor.

        Parameters
        ----------
        sim : Simulation instance
            The simulation object to which this Cantor catalogue belongs.
        index : int, optional
            The output index of this catalogue. If None (default), the
            object refers to global quantities.
        snap_type : str, optional
            Whether this Cantor catalogue belongs to a snapshot
            ('snap', default), or snipshot ('snip').
        cantor_dir : str, optional
            The name of the Cantor catalogue directory. If None (default),
            it is loaded from the Simulation instance sim.
        read_range : (int, int) or None, optional
            Read only elements from the first up to *but excluding* the
            second entry in the tuple (in dimension 0). If None (default),
            read the entire file. Ignored if read_index is provided.
        read_index : int or np.array(int) or None, optional
            Read only the elements in read_index (in dimension 0). If int,
            a single element is read, and the first dimension truncated. If
            an array is provided, the elements between the lowest and
            highest index are read and the output then masked to the
            exact elements. If None (default), everything is read.
        """
        # Set up paths to catalogue files
        if cantor_dir is None:
            cantor_dir = sim.cantor_dir
        if index is not None:
            self.catloc = cantor_dir + f'Cantor_{index:03d}.hdf5'
            if not os.path.isfile(self.catloc):
                print(f"Cantor catalogue {index} not found!")
                set_trace()
                self.catloc = None
                return
            self.idloc = cantor_dir + f'Cantor_{index:03d}_IDs.hdf5'
            if not os.path.isfile(self.idloc):
                self.idloc = None
            self.radloc = cantor_dir + f'Cantor_{index:03d}_Radii.hdf5'
            if not os.path.isfile(self.radloc):
                self.radloc = None
        else:
            print("Setting up cross-snapshot functionality is still "
                  "in progress. Please come back later...")
            set_trace()

        self.group = group
        self.read_index = read_index
        if read_index is None:
            self.read_range = read_range
        elif isinstance(read_index, int):
            self.read_range = (read_index, read_index+1)
        else:
            self.read_range = (np.min(read_index), np.max(read_index)+1)

    def __getattr__(self, name):
        """Provide a way to handle requests for new attributes."""
        if name == 'IDs':
            if self.idloc is None:
                print("IDs file not found -- cannot load!")
                set_trace()
            setattr(self, name, hy.hdf5.read_data(
                self.idloc, name, self.read_range, self.read_index))
        elif name == 'Radii':
            if self.radloc is None:
                print("Radius file not found -- cannot load!")
                set_trace()
            setattr(self, name, hy.hdf5.read_data(
                self.radloc, name, self.read_range, self.read_index))

        # 'Normal' situation: read from actual catalogue.
        else:
            # Exception: reading extra catalogue entries
            if (name.startswith('Extra') and name != 'Extra__ExtraIDs'
                and self.read_index is not None):

                # Read relevant extra entries
                read_index = self.Extra__ExtraIDs
                matched = np.nonzero(read_index >= 0)[0]
                if isinstance(self.read_index, int):
                    read_index = np.array([read_index])
                data = hy.hdf5.read_data(
                    self.catloc, self.group + '/' + name.replace('__', '/'),
                    read_index=read_index[matched])
                if data is None:
                    print(f"Could not find requested data '{name}'")
                    setattr(self, name, data)
                    return None

                # Re-align the data to main catalogue
                data_shape = list(data.shape)
                if isinstance(self.read_index, int):
                    data_shape[0] = 1
                else:
                    data_shape[0] = len(self.read_index)
                data_out = np.zeros(data_shape, data.dtype) - 1
                data_out[matched, ...] = data
                data = data_out
                if isinstance(self.read_index, int):
                    if data.ndim == 1:
                        data = data[0]
                    else:
                        data = data[0, ...]
                setattr(self, name, data)

            # Normal situation: reading main catalogue entries
            else:
                data = hy.hdf5.read_data(
                self.catloc, self.group + '/' + name.replace('__', '/'),
                self.read_range, self.read_index)

                if data is None:
                    print(f"Could not find requested data '{name}'")
                setattr(self, name, data)

        return getattr(self, name)


