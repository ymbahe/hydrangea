"""Provides SplitFile class to read distributed HDF5 data sets."""

import os
import numpy as np
import h5py as h5

from pdb import set_trace
import hydrangea.hdf5 as hd
import hydrangea.tools as ht
from hydrangea.reader_base import ReaderBase
import time


class SplitFile(ReaderBase):
    """Class to read data sets that are split over multiple files.

    Attributes
    ----------
    file_name : str
        Name of one of the files in this split collection.
    group_name : str or None
        Base group to read data from (None if not specified)
    num_elem : int or None
        Number of elements in the selected group category (None if it
        could not be determined or is not applicable).
    num_files : int or None
        Number of files in the collection (None if not determined).
    """

    def __init__(self, file_name, group_name=None, part_type=None,
                 sim_type='Eagle', verbose=False):
        """Initialize a file collection to read from.

        Parameters
        ----------
        file_name : str
            Any one of the files in the collection to read.
        group_name : str, optional
            Name of base group to read from. If not provided, it is
            assumed to be 'PartType[x]', where x is the particle type
            (which must then be supplied).
        part_type : int, optional
            The (numerical) particle type; only relevant for sn[a/i]pshots.
            If not provided, it is inferred from group_name, which must
            then be supplied.
        sim_type : str, optional
            Type of simulation to which this data belongs:
            'Eagle' (default) or 'Illustris'.
        verbose : bool, optional
            Provide extended output (default: False).
        """
        self.verbose = verbose

        # First: check if the file exists...!
        if not os.path.isfile(file_name):
            print("-----------------------------------------------------")
            print("Error [" + file_name + "]")
            print("This file does not exist. Please supply correct name!")
            print("-----------------------------------------------------")
            return
        self.file_name = file_name
        self.sim_type = sim_type

        # Set up group name, if needed from particle type
        if part_type is not None:
            self._decode_ptype(part_type)
        else:
            self.group_name = group_name

        if verbose:
            print("Prepared reading from '{:s}'..."
                  .format(self.group_name.upper()))

        # For now, don't build offset list
        self.file_offsets = None

    def read_data(self, dataset_name, verbose=False, astro=True,
                  return_conv=False):
        """Read a specified data set from the file collection.

        Parameters
        ----------
        dataset_name : str
            The name of the data set to read, including possibly containing
            groups, but *not* the main group specified in the instantiation.
        astro : bool
            Attempt conversion to `astronomical' units where necessary
            (pMpc, 10^10 M_sun, km/s). This is ignored for dimensionless
            quantities. Default: True
        verbose : bool
            Provide more or less useful messages during reading.
            Default: False

        Returns
        -------
        data : np.array
            Array containing the specified data.
        astro_conv : float or None
            Conversion factor to astronomical units; returned only if
            return_conv is True.

        """
        # Check that setup has been done properly
        if self.num_files is None:
            print("I don't know how many files to read...")
            set_trace()

        # Read individual files
        if self.num_elem is not None:
            if self.file_offsets is not None:
                data_out = self._read_files_direct(dataset_name)
            else:
                data_out = self._read_files_sequentially(dataset_name)
        else:
            # Use HAC if we don't know how many elements to read in total
            data_out = self._read_files_hac(dataset_name)

        # Apply corrections if needed
        if astro or return_conv:
            astro_conv = self.get_astro_conv(dataset_name)
            if astro and astro_conv is not None:
                data_out *= astro_conv

        if return_conv:
            return data_out, astro_conv
        else:
            return data_out

    def _read_files_direct(self, dataset_name, verbose=True):
        """Read data set from files using pre-established offset list."""
        data_out = None
        for ifile in range(self.num_files):
            num_exp = self.file_offsets[ifile+1] - self.file_offsets[ifile]
            if num_exp > 0:
                num_read, data_out = self._read_file(
                    ifile, dataset_name, data_out,
                    offset=self.file_offsets[ifile], verbose=verbose)
                if num_read != num_exp:
                    print("Read {:d} elements, but expected {:d}!"
                          .format(num_read, num_exp))
        print("")
        return data_out

    def _read_files_sequentially(self, dataset_name, verbose=True):
        """Read data from files in sequential order."""
        data_out = None
        offset = 0

        # Read data into data_out, which is cycled through _read_file()
        for ifile in range(self.num_files):
            num_read, data_out = self._read_file(
                ifile, dataset_name, data_out, offset=offset, verbose=verbose)
            offset += num_read

        # Make sure we read right total number
        if offset != self.num_elem:
            print("Read wrong number of elements for '{:s}'."
                  .format(dataset_name))
            set_trace()

        print("")   # Ends no-newline sequence inside _read_file()
        return data_out

    def _read_files_hac(self, dataset_name, threshold=int(1e9)):
        """Read data set from files into output, using HAC."""
        data_out = None
        data_part = None

        for ifile in range(self.num_files):
            num, data_part = self._read_file(ifile, dataset_name, None,
                                             self_only=True)

            # Initialize output and stack, or append part to stack
            if data_out is None:
                data_out = np.copy(data_part)
                data_stack = np.copy(data_part)
            else:
                data_stack = np.concatenate((data_stack, data_part))

            # Combine to full list if critical size reached:
            if len(data_stack) > threshold:
                print("Update full list...")
                data_out = np.concatenate((data_out, data_stack))
                empty_shape = list(data_stack.shape)
                empty_shape[0] = 0
                data_stack = np.empty(empty_shape, data_stack.dtype)

        print("")  # Ends no-newline sequence inside _read_file()

        # Need to do final concatenation after loop ends
        if data_stack:
            data_out = np.concatenate((data_out, data_stack))

        return data_out

    def _read_file(self, ifile, dataset_name, data_out, offset=0,
                   verbose=True, read_range=None, num_out=None,
                   self_only=False):
        """Read specified data set from one file into output.

        This is the low-level reading routine called by the three
        'driver' functions _read_files_... above.

        Parameters
        ----------
        ifile : int
            The index of the file to read.
        dataset_name : str
            The name of the data set to read, including possibly containing
            groups but *not* the main base group (e.g. 'PartType0').
        data_out : np.array or None
            The output array. If None, it is initialized internally.
        offset : int, optional
            The offset into the output array to write the data to.
            Default is 0, i.e. fill output array from its beginning.
        verbose : bool, optional
            Print progress messages.
        read_range : int [start, end] or None, optional
            The range of the data to read. If None (default), read all.
            THIS IS STILL EXPERIMENTAL AND NOT PROPERLY IMPLEMENTED.
        num_out : int or None, optional
            The number of elements to allocate in the output array, if
            it is initialized internally. If None (default), use self.numElem.
        self_only : bool, optional
            If output is internally initialized, only allocate enough
            elements for data from this file (default: False).

        Returns
        -------
        length : int
            The number of elements that were read.
        data_out : np.array or None
            The updated (or newly initialized) output array.
        """

        if verbose:
            print(str(ifile) + " ", end="", flush=True)

        # Form current file name and check it exists
        file_name = self._swap_file_name(self.file_name, ifile)
        if not os.path.isfile(file_name):
            print("\nError: did not find expected file {:d}."
                  .format(ifile), flush=True)
            set_trace()

        f = h5.File(file_name, 'r')
        full_dataset_name = self.group_name + '/' + dataset_name

        # Not all files contain all groups, so need to check explicitly
        try:
            dataSet = f[full_dataset_name]
        except KeyError:
            if verbose:
                print("No data found on file {:d}!" .format(ifile))
            return 0, data_out

        # With read range, length is pre-determined
        if read_range is not None:
            length = read_range[1] - read_range[0]

        if data_out is None:
            if offset > 0:
                print("Warning: initializing output although offset is {:d}."
                      .format(offset))
            shape = list(dataSet.shape)   # Need to modify, so must be list

            if read_range is None:
                length = shape[0]

            if not self_only:
                if num_out is None:
                    shape[0] = self.num_elem
                else:
                    shape[0] = num_out
            data_out = np.empty(shape, dataSet.dtype)
            dataSet.read_direct(data_out,
                                dest_sel=np.s_[offset:offset+length])
        else:
            if read_range is None:
                length = dataSet.len()
            dataSet.read_direct(data_out,
                                dest_sel=np.s_[offset:offset+length])

        return length, data_out

    def _decode_ptype(self, ptype):
        """Identify supplied particle type"""

        if isinstance(ptype, int):
            self.group_name = 'PartType{:d}' .format(ptype)
        elif isinstance(ptype, str):
            if ptype.upper() == 'GAS':
                self.group_name = 'PartType0'
            elif ptype.upper() == 'DM':
                self.group_name = 'PartType1'
            elif ptype.upper() in ['STARS', 'STAR']:
                self.group_name = 'PartType4'
            elif ptype.upper() in ['BH', 'BHS', 'BLACKHOLE', 'BLACKHOLES',
                                   'BLACK_HOLES', 'BLACK_HOLE']:
                self.group_name = 'PartType5'
            else:
                print("Unrecognized particle type name '{:s}'"
                      .format(ptype))
                set_trace()
        else:
            print("Unrecognized particle type format")
            set_trace()

    @property
    def num_elem(self):
        """Find out how many output elements there are in total."""
        if '_num_elem' not in dir(self):
            self._num_elem = self._count_elements()
        return self._num_elem

    def _count_elements(self):
        """Count number of elements to read."""
        num_elem = None   # Placeholder for "not known"
        if self.group_name is None:
            return

        # Break file name into base and running sequence number
        real_file_name = os.path.split(self.file_name)[1]
        file_name_parts = real_file_name.split('_')

        if self.group_name.startswith('PartType'):
            pt_index = int(self.group_name[8])
        else:
            pt_index = None   # For checking

        # Deal with particle catalogue files
        if (file_name_parts[0] in ['snap', 'snip', 'partMags',
                                   'eagle_subfind_particles']):
            if self.verbose:
                print("Particle catalogue detected!")
            # Need to extract particle index to count (mag --> stars!)
            if self.group_name.startswith('PartType'):
                num_elem = self._count_elements_snap(pt_index)
            elif file_name_parts[0] == 'partMags':
                num_elem = self._count_elements_snap(4)
            else:
                return    # Can't determine element numbers then

        # Deal with subfind catalogue files
        elif len(file_name_parts) >= 2:
            if "_".join(file_name_parts[:3]) == 'eagle_subfind_tab':
                if self.verbose:
                    print("Subfind catalogue detected!")
                if self.group_name == 'FOF':
                    num_elem = self._count_elements_sf_fof()
                elif self.group_name == 'Subhalo':
                    num_elem = self._count_elements_sf_subhalo()
                elif self.group_name == 'IDs':
                    num_elem = self._count_elements_sf_ids()
                else:
                    return

        # In all other cases, don't know how to pre-determine
        # element count, so need to find out from individual files (later).
        else:
            return

        # Final piece: deal with possibility of int32 overflow in numbers
        # (it does happen somewhere...)
        if num_elem < 0:
            num_elem += 4294967296
        if self.verbose:
            print("   ... (detected {:d} elements) ..."
                  .format(num_elem), flush=True)
        return num_elem

    def _count_elements_snap(self, pt_index):
        """Count number of particles of specified (int) type."""
        return hd.read_attribute(
            self.file_name, 'Header', 'NumPart_Total', require=True)[pt_index]

    def _count_elements_sf_fof(self):
        """Count number of FOF groups in Subfind catalogue."""
        return hd.read_attribute(
            self.file_name, 'Header', 'TotNgroups', require=True)

    def _count_elements_sf_subhalo(self):
        """Count number of subhaloes in Subfind catalogue."""
        return hd.read_attribute(
            self.file_name, 'Header', 'TotNsubgroups', require=True)

    def _count_elements_sf_ids(self):
        """Count number of particle IDs in Subfind catalogue."""
        return hd.read_attribute(
            self.file_name, 'Header', 'TotNids', require=True)

    @property
    def num_files(self):
        """Count number of files to read."""
        if '_num_files' not in dir(self):
            if self.sim_type == 'Eagle':
                self._num_files = hd.read_attribute(
                    self.file_name, 'Header', 'NumFilesPerSnapshot',
                    require=True)
            elif self.sim_type == 'Illustris':
                self._num_files = hd.read_attribute(
                    self.file_name, 'Header', 'NumFiles', require=True)
            else:
                print("Unknown simulation type '{:s}'." .format(self.sim_type))
                set_trace()
        return self._num_files