"""Provides a base class for specialized reader classes."""

import numpy as np
import h5py as h5
import hydrangea.tools as ht
import hydrangea.hdf5 as hd
import hydrangea.crossref as hx
from hydrangea.split_file import SplitFile
from pdb import set_trace  # [Will be used later]


class ReaderBase:
    """Base for reading classes (SplitFile and ReadRegion).

    It contains common functionality, in particular conversion to
    astronomical units, on-the-fly luminosity calculation, and
    cross-matching of particle IDs to subhalo catalogues.
    """

    def get_astro_conv(self, dataset_name):
        """Get the conversion factor to astronomical units for a data set.

        Parameters
        ----------
        dataset_name : str
            The data set for which to compute the conversion factor.

        Returns
        -------
        conv : float or None
            Conversion factor (None if it could not be determined).
        """
        try:
            f = h5.File(self.file_name, 'r')
            dSet = f[self.base_group + '/' + dataset_name]

            hscale_exponent = dSet.attrs["h-scale-exponent"]
            ascale_exponent = dSet.attrs["aexp-scale-exponent"]

            header = f["/Header"]
            aexp = header.attrs["ExpansionFactor"]
            h_hubble = header.attrs["HubbleParam"]
            f.close()
        except(KeyError):
            return None

        conv_astro = aexp**ascale_exponent * h_hubble**hscale_exponent
        return conv_astro

    @property
    def aexp(self):
        """Get the expansion factor of the data set."""
        if '_aexp' not in dir(self):
            self._aexp = hd.read_attribute(self.file_name, 'Header',
                                           'ExpansionFactor', require=True)
        return self._aexp

    @property
    def redshift(self):
        """Get the redshift of the data set."""
        if '_zred' not in dir(self):
            self._zred = ht.aexp_to_time(self.aexp, 'zred')
        return self._zred

    @property
    def time(self):
        """Get the age of the Universe of the data set [Gyr]."""
        if '_age' not in dir(self):
            self._time = ht.aexp_to_time(self.aexp, 'age')
        return self._time

    @property
    def lookback_time(self):
        """Get the lookback time to the data set from z = 0 [Gyr]."""
        if '_lbt' not in dir(self):
            self._lbt = ht.aexp_to_time(self.aexp, 'lbt')
        return self._lbt

    @property
    def m_dm(self):
        """Get the DM particle mass of a snapshot."""
        if '_m_dm' not in dir(self):
            if not self.base_group.startswith('PartType'):
                print("*** Looks like you are trying to load m_dm for "
                      "non-snapshot! ***")
                set_trace()
            self._m_dm = ht.get_m_dm(self.file_name, astro=self.astro)
        return self._m_dm

    @property
    def m_baryon(self):
        """Get the initial baryon mass of a snapshot."""
        if '_m_baryon' not in dir(self):
            if not self.base_group.startswith('PartType'):
                print("*** Looks like you are trying to load m_baryon for "
                      "non-snapshot! ***")
                set_trace()
            self.m_baryon = ht.get_m_baryon(self.file_name, astro=self.astro)
        return self._m_baryon


    def get_time(self, timeType='aexp'):
        """Retrieve various possible time stamps of the output.

        Parameters
        ----------
        timeType : string
            Specified the 'time' flavour to be returned. Options are:
                'aexp' [default]: Expansion factor
                'zred': Redshift
                'age': Age of the Universe [Gyr]
                'lbt': Lookback time from z = 0 [Gyr]

        Returns
        -------
        time : float
            The requested time stamp.
        """
        return ht.aexp_to_time(self.aexp, timeType)

    def _swap_file_name(self, base_name, number):
        """Form sequentially numbered file names (convenience function)."""
        name_parts = base_name.split('.')
        name_parts[-2] = str(number)
        new_name = '.'.join(name_parts)
        return new_name

    def __getattr__(self, name):
        """Provide a way to handle requests for non-existing attributes."""
        self._print(2, "Requested attribute '{:s}' does not yet exist..."
                    .format(name))
        if '__' in name:
            name_actual = name.replace('__', '/')
        else:
            name_actual = name
        if name_actual == 'SubhaloIndex':
            if 'SubhaloIndex' not in dir(self):
                self.SubhaloIndex = self.find_group(group_type='subfind')
            return self.SubhaloIndex
        elif name_actual == 'GroupIndex':
            if 'GroupIndex' not in dir(self):
                self.GroupIndex = self.find_group(group_type='fof')
            return self.GroupIndex
        else:
            return self.read_data(name_actual, store=None, trial=True)

    def _print(self, threshold, *args, **kwargs):

        if isinstance(threshold, int):
            verbose = self.verbose
            thresh = threshold
        else:
            thresh = threshold[0]
            if len(threshold) == 1 or threshold[1] is None:
                verbose = self.verbose
            else:
                verbose = threshold[1]

        if verbose >= thresh:
            print(*args, **kwargs)

    def find_group(self, ids=None, group_type="subfind",
                   return_matched=False):
        """On-the-fly retrieval of particle group indices.

        This is achieved by cross-matching particle IDs between the
        snapshot catalogue and the corresponding Subfind catalogue.

        Parameters
        ----------
        ids : np.array(int), optional
            The particle IDs to match. If None (default), match all
            particles of this data set.
        group_type : str, optional
            Type of group to find membership for. Possible values are
            'subfind' (default) or 'fof'.
        return_matched : bool, optional
            Also find and return elements that are actually in a target
            group block (default: False).

        Returns
        -------
        shi : np.array(int)
            The subhalo index for each particle (-1 for particles which are
            not in a subhalo).
        in_sh : np.array(int)
            The indices of particles that could be located in a subhalo
            (i.e. those whose shi is >= 0).

        Note
        ----
        This is a convenience function to emulate a (non-existing)
        'subhalo index' data set in snapshot catalogues. Depending on
        circumstances, other approaches may be faster.
        """
        if ids is None:
            ids = self.ParticleIDs

        if group_type.lower() == 'subfind':
            group_ids = SplitFile(self.subfind_file, 'IDs')
            group_cat = SplitFile(self.subfind_file, 'Subhalo')
            index_in_ids = hx.find_id_indices(ids, group_ids.ParticleID)[0]
            offset = group_cat.SubOffset
            group = ht.ind_to_block(index_in_ids, offset, group_cat.SubLength)
        elif group_type.lower() == 'fof':
            group_ids = SplitFile(self.subfind_file, 'IDs')
            group_cat = SplitFile(self.subfind_file, 'FOF')
            index_in_ids = hx.find_id_indices(ids, group_ids.ParticleID)[0]
            offset = group_cat.GroupOffset
            group = ht.ind_to_block(
                index_in_ids, offset, group_cat.GroupLength)
        else:
            print(f"Invalid group_type ('{group_type}')")
            set_trace()

        if return_matched:
            ind_matched = np.nonzero(
                (group >= 0) & (group < len(offset)))[0]
            return group, ind_matched
        else:
            return group

    @property
    def subfind_file(self):
        """Subfind catalogue file belonging to a snapshot."""
        if '_subfind_file' not in dir(self):
            print("Need to set a corresponding subfind file!")
            set_trace()
        return self._subfind_file

    @subfind_file.setter
    def subfind_file(self, cat_file):
        """Set a subfind catalogue file to associate with this file."""
        self._subfind_file = cat_file

    @property
    def cantor_file(self):
        """Cantor catalogue file belonging to a snapshot."""
        if '_cantor_file' not in dir(self):
            print("Need to set a corresponding cantor file first!")
            set_trace()
        return self._cantor_file

    @cantor_file.setter
    def cantor_file(self, cat_file):
        """Set a cantor catalogue file to associate with this file."""
        self._subfind_file = cat_file
