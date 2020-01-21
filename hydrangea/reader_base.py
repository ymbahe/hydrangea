"""Provides a base class for specialized reader classes."""

import numpy as np
import h5py as h5
import hydrangea.tools as ht
import hydrangea.units as hu
import hydrangea.hdf5 as hd

from pdb import set_trace  # [Will be used later]


class ReaderBase:
    """Base for reading classes (SplitFile and ReadRegion).

    It contains common functionality, in particular conversion to
    astronomical units, on-the-fly luminosity calculation, and
    cross-matching of particle IDs to subhalo catalogues.
    """

    def get_unit_conversion(self, dataset_name, units_name):
        """Get appropriate factor to convert data units to other system."""
        if units_name == 'data':
            data_to_other = 1
        elif units_name == 'clean':
            data_to_other = self.get_clean_factor(self, dataset_name)
        elif units_name == 'astro':
            data_to_other = (self.get_clean_factor(self, dataset_name)
                             * self.get_astro_factor(self, dataset_name))
        elif units_name == 'cgs':
            data_to_other = (self.get_clean_factor(self, dataset_name)
                             * self.get_cgs_factor(self, dataset_name))
        elif units_name == 'si':
            data_to_other = (self.get_clean_factor(self, dataset_name)
                             * self.get_si_factor(self, dataset_name))
        else:
            print(f"Unknown unit system requested: '{units_name}'!")
            set_trace()

        return data_to_other

    def get_clean_factor(self, dataset_name):
        """Get the conversion factor to 'clean data units' (no h or a)."""
        f = h5.File(self.file_name, 'r')
        dSet = f[self.base_group + '/' + dataset_name]

        try:
            hscale_exponent = dSet.attrs["h-scale-exponent"]
            ascale_exponent = dSet.attrs["aexp-scale-exponent"]

            header = f["/Header"]
            aexp = header.attrs["ExpansionFactor"]
            h_hubble = header.attrs["HubbleParam"]
        except(KeyError):
            f.close()
            return None
        return aexp**ascale_exponent * h_hubble**hscale_exponent

    def get_astro_factor(self, dataset_name):
        """Get the conversion factor to 'astronomical units'."""
        data_to_cgs_factor = self.get_cgs_factor(self, dataset_name)
        if data_to_cgs_factor is None:
            return None

        dimensions = hu.get_dimensions(self.base_group, dataset_name)
        cgs_to_astro_factor = (hu.si_to_astro_factor(dimensions)
                               / hu.si_to_cgs_factor(dimensions))
        return data_to_cgs_factor * cgs_to_astro_factor

    def get_cgs_factor(self, dataset_name):
        """Get the conversion factor to CGS units (boo, hiss)."""
        return hd.read_attribute(self.file_name, dataset_name, 'CGS_Factor')

    def get_si_factor(self, dataset_name):
        """Get the conversion factor to SI units."""
        data_to_cgs_factor = self.get_cgs_factor(self, dataset_name)
        if data_to_cgs_factor is None:
            return None

        dimensions = hu.get_dimensions(self.base_group, dataset_name)
        cgs_to_si_factor = 1 / hu.si_to_cgs_factor(dimensions)
        return data_to_cgs_factor * cgs_to_si_factor

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
            self._m_dm = ht.get_m_dm(self.file_name, units=self.units)
        return self._m_dm

    @property
    def m_baryon(self):
        """Get the initial baryon mass of a snapshot."""
        if '_m_baryon' not in dir(self):
            if not self.base_group.startswith('PartType'):
                print("*** Looks like you are trying to load m_baryon for "
                      "non-snapshot! ***")
                set_trace()
            self.m_baryon = ht.get_m_baryon(self.file_name, units=self.units)
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
