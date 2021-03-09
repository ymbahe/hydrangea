"""Provides a base class for specialized reader classes."""

import h5py as h5
import hydrangea.tools as ht
import hydrangea.units as hu
import hydrangea.hdf5 as hd
import numpy as np

from pdb import set_trace  # [Will be used later]


class ReaderBase:
    """Base for reading classes (SplitFile and ReadRegion).

    It contains common functionality, in particular conversion to
    specified unit systems and cross-matching of particle IDs to
    subhalo catalogues.
    """

    def get_unit_conversion(self, dataset_name, units_name):
        """Get appropriate factor to convert data units to other system.

        Parameters
        ----------
        dataset_name : str
            The name of the data set for which to obtain the conversion
            factor (including possible containing groups, but not the
            base group).
        units_name : str
            The unit system to calculate the conversion factor for.
            Options are (case-insensitive):

            - ``'data'`` --> Exactly as stored in file (i.e. no conversion)
            - ``'clean'`` --> As in file, but without `a` and `h` factors
            - ``'astro'`` --> Astronomically useful units (e.g. M_Sun, pMpc)
            - ``'si'`` --> SI units
            - ``'cgs'`` --> CGS units

        Returns
        -------
        data_to_other : float
            The conversion factor for the specified unit system.
            The 'raw' data as read from the file(s) must be multiplied
            with this value to obtain the magnitude in the target system.

        Note
        ----
        In particular in SI and CGS units, overflow issues may occur for
        32-bit floats, because these have a maximum value of ~1e39.

        Examples
        --------
        >>> import hydrangea as hy
        >>> snap_file = hy.objects.Simulation(index=0).get_snap_file(29)
        >>> stars = hy.SplitFile(snap_file, part_type=4)  # or hy.ReadRegion
        >>> stars.get_unit_conversion('Mass', 'astro')
        14755791648.22193

        >>> stars.get_unit_conversion('CentreOfPotential', 'SI')
        4.553162166150214e+22
        """
        if units_name.lower() == 'data':
            data_to_other = 1
        elif units_name.lower() == 'clean':
            data_to_other = self.get_data_to_clean_factor(dataset_name)
        elif units_name.lower() == 'astro':
            data_to_other = self.get_data_to_astro_factor(dataset_name)
        elif units_name.lower() == 'cgs':
            data_to_other = self.get_data_to_cgs_factor(dataset_name)
        elif units_name.lower() == 'si':
            data_to_other = self.get_data_to_si_factor(dataset_name)
        else:
            print(f"Unknown unit system requested: '{units_name}'!")
            set_trace()

        return data_to_other

    def get_data_to_clean_factor(self, dataset_name):
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

    def get_data_to_cgs_factor(self, dataset_name):
        """Get the conversion factor to CGS units (boo, hiss)."""
        data_to_clean_factor = self.get_data_to_clean_factor(dataset_name)
        if data_to_clean_factor is None:
            return None

        return (hd.read_attribute(self.file_name,
                                  self.base_group + '/' + dataset_name,
                                  'CGSConversionFactor')
                * data_to_clean_factor)

    def get_data_to_astro_factor(self, dataset_name):
        """Get conversion factor from data to astro."""
        data_to_cgs_factor = self.get_data_to_cgs_factor(dataset_name)
        if data_to_cgs_factor is None:
            return None

        cgs_to_astro_factor = self.get_cgs_to_astro_factor(dataset_name)
        if cgs_to_astro_factor is None:
            return None
        else:
            return data_to_cgs_factor * cgs_to_astro_factor

    def get_data_to_si_factor(self, dataset_name):
        """Get conversion factor from data to SI."""
        data_to_cgs_factor = self.get_data_to_cgs_factor(dataset_name)
        if data_to_cgs_factor is None:
            return None

        cgs_to_si_factor = self.get_cgs_to_si_factor(dataset_name)
        if cgs_to_si_factor is None:
            return None
        else:
            return data_to_cgs_factor * cgs_to_si_factor

    def get_cgs_to_astro_factor(self, dataset_name):
        """Get the conversion factor from CGS to 'astronomical units'."""
        dimensions = hu.get_dimensions(self.base_group, dataset_name)
        return hu.get_unit_conversion_factor(dimensions, 'cgs', 'astro')

    def get_cgs_to_si_factor(self, dataset_name):
        """Get the conversion factor from CGS to SI units."""
        dimensions = hu.get_dimensions(self.base_group, dataset_name)
        return hu.get_unit_conversion_factor(dimensions, 'cgs', 'si')

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
        """Expansion factor of the data set."""
        if '_aexp' not in dir(self):
            self._aexp = hd.read_attribute(self.file_name, 'Header',
                                           'ExpansionFactor', require=True)
        return self._aexp

    @property
    def redshift(self):
        """Redshift of the data set."""
        if '_zred' not in dir(self):
            self._zred = ht.aexp_to_time(self.aexp, 'zred')
        return self._zred

    @property
    def time(self):
        """Age of the Universe of the data set [Gyr]."""
        if '_age' not in dir(self):
            self._time = ht.aexp_to_time(self.aexp, 'age')
        return self._time

    @property
    def lookback_time(self):
        """Lookback time to the data set from z = 0 [Gyr]."""
        if '_lbt' not in dir(self):
            self._lbt = ht.aexp_to_time(self.aexp, 'lbt')
        return self._lbt

    @property
    def m_dm(self):
        """DM particle mass (only for snap-/snipshots)."""
        if '_m_dm' not in dir(self):
            if not self.base_group.startswith('PartType'):
                print("*** Looks like you are trying to load m_dm for "
                      "non-snapshot! ***")
                set_trace()
            self._m_dm = ht.get_m_dm(self.file_name, units=self.read_units)
        return self._m_dm

    @property
    def m_baryon(self):
        """Initial baryon mass (only for snap-/snipshots)."""
        if '_m_baryon' not in dir(self):
            if not self.base_group.startswith('PartType'):
                print("*** Looks like you are trying to load m_baryon for "
                      "non-snapshot! ***")
                set_trace()
            self.m_baryon = ht.get_m_baryon(self.file_name, 
                                            units=self.read_units)
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
        elif name_actual == 'Mass' and self.base_group == 'PartType1':
            if 'Mass' not in dir(self):
                self.Mass = np.zeros(self.num_elem) + self.m_dm
            return self.Mass
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
        """Subfind catalogue file associated to a snapshot.

        This must be set explicitly by the user. It can point to the
        output's own subfind file (if it exists), but does not need to:
        in the latter case, it allows matching particles to structures
        at another point in time.
        """
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
