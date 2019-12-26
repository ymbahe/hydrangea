import numpy as np
import h5py as h5
import hydrangea.crossref as hx

class ReaderBase:
    """Base for reading classes (SplitFile and ReadRegion).
    
    It contains common functionality, in particular conversion to 
    astronomical units, on-the-fly luminosity calculation, and 
    cross-matching of particle IDs to subhalo catalogues.
    """

    def get_astro_conv(self, dataSetName):
        """Get the conversion factor to astronomical units for a data set.

        Parameters
        ----------
        dataSetName : str
            The data set for which to compute the conversion factor.

        Returns
        -------
        conv : float or None
            Conversion factor (None if it could not be determined).
        """

        try:
            f = h5.File(self.fileName, 'r')
            dSet = f[self.baseGroup + '/' + dataSetName]
            
            hscale_exponent = dSet.attrs["h-scale-exponent"]
            ascale_exponent = dSet.attrs["aexp-scale-exponent"]
    
            header = f["/Header"]
            aexp = header.attrs["ExpansionFactor"]
            h_hubble = header.attrs["HubbleParam"]
            f.close()

        except:
            return None
    
        conv_astro = aexp**ascale_exponent * h_hubble**hscale_exponent
        return conv_astro

    def _swap_file_name(self, baseName, number):
        """Convenience function to update base filename"""

        nameParts = baseName.split('.')
        nameParts[-2] = str(number)
        resName = '.'.join(nameParts)
        return resName
    
    def find_subhalo(self, ids=None, return_matched=False):
        """On-the-fly retrieval of particle subhalo indices, found by
        cross-matching particle IDs.
        
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
        'subhalo index' data set in sn(a/i)pshot catalogues. Depending on
        circumstances, other approaches may be faster.
        """

        pass
        
