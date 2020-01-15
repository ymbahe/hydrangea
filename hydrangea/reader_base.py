"""Provides a base class for specialized reader classes."""

import h5py as h5
import hydrangea.tools as ht
import hydrangea.hdf5 as hd
import hydrangea.crossref as hx  # [Will be used later]
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
            return self.read_data(name_actual, store=None, trial=True)
        return self.read_data(name, store=None, trial=True)

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

    def find_subhalo(self, ids=None, return_matched=False):
        """On-the-fly retrieval of particle subhalo indices.

        This is achieved by cross-matching particle IDs between the
        snapshot catalogue and the corresponding Subfind catalogue.

        Parameters
        ----------
        ids : np.array(int), optional
            The previously read particle IDs. If None (default), these will
            be read in internally, at additional computational cost.

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

        This functionality is not yet implemented.
        """
        pass