# Blank file

import os
from pdb import set_trace
import hydrangea.hdf5 as hd
import hydrangea.tools as ht

class SplitFile:
    """Class to read data sets that are split over multiple files.

    Attributes
    ----------
    fileName : str
        Name of one of the files in this split collection.
    groupName : str or None
        Base group to read data from (None if not specified)
    numElem : int or None
        Number of elements in the selected group category (None if it 
        could not be determined or is not applicable).
    numFiles : int or None
        Number of files in the collection (None if not determined).
    """
        
    def __init__(self, fileName, groupName=None, ptype=None, verbose=False,
                 simType='Eagle'):
                
        # First: check if the file exists...!
        if not os.path.isfile(fileName):
            print("-----------------------------------------------------")
            print("Error [" + fileName + "]")
            print("This file does not exist. Please supply correct name!")
            print("-----------------------------------------------------")
            return 
        self.fileName = fileName
        
        # Set up group name, if needed from particle type
        if ptype is not None:
            self._decode_ptype(ptype)
        else:
            self.groupName = groupName
  
        # Find the total number of elements in the output array
        self._count_elements()         # Stored in self.numElem
        self._count_files(simType)     # Stored in self.numFiles
                
        if verbose:
            print("Prepared reading from '{:s}'..."
                  .format(self.groupName.upper()))

    def read(self, dataSetName, verbose=False, astro=True,
                 return_conv=False):
        """Read a specified data set from the file collection.
        
        Parameters
        ----------
        dataSetName : str
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
        if self.numFiles is None:
            print("I don't know how many files to read...")
            set_trace()
        
        # Set up output array if needed
        if setup_before_read:
            data_out = self._setup_output(dataSetName)
        else:
            data_out = None

        # Read individual files
        if self.numElem is not None:
            for ifile in range(self.numFiles):
                self._read_file(ifile, dataSetName, data_out)
        else:
            # Use HAC if we don't know how many elements to read in total
            self._read_files_hac(dataSetName, data_out)
                
        # Apply corrections if needed
        if astro or return_conv:
            astro_conv = self.get_astro_conv(dataSetName)
            if astro and astro_conv is not None:
                data_out *= astro_conv

        if return_conv:
            return data_out, astro_conv
        else:
            return data_out
    
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
        aExp = hd.read_attribute(
            self.fileName, 'Header', 'ExpansionFactor', require=True)
        return ht.aexp_to_time(timeType)

    def get_astro_conv(self, dataSetName):
        """Retrieve the conversion factor to 'astronomical' units 
        for a specified data set.

        Note
        ----
        Returns None if this is an inconvertible quantity (e.g. IDs).
        """
        return ht.get_astro_conv(self.fileName, dataSetName)

    def _setup_output(self, dataSetName):
        """Set up the output array for specified data set."""
        pass

    def _read_file(self, ifile, dataSetName, data_out):
        """Read a specified data set from a specified file into output."""
        pass

    def _read_files_hac(self, dataSetName, data_out):
        """Read data set from files into output, using HAC."""
        pass
        
    def _decode_ptype(self, ptype):
        """Identify supplied particle type"""
        
        if isinstance(ptype, int):
            self.groupName = 'PartType{:d}' .format(ptype)
        elif isinstance(ptype, str):
            if ptype.upper() == 'GAS':
                self.groupName = 'PartType0'
            elif ptype.upper() == 'DM':
                self.groupName = 'PartType1'
            elif ptype.upper() in ['STARS', 'STAR']:
                self.groupName = 'PartType4'
            elif ptype.upper() in ['BH', 'BHS', 'BLACKHOLE', 'BLACKHOLES',
                                   'BLACK_HOLES', 'BLACK_HOLE']:
                self.groupName = 'PartType5'
            else:
                print("Unrecognized particle type name '{:s}'"
                      .format(ptype))
                set_trace()
        else:
            print("Unrecognized particle type format")
            set_trace()

    def _count_elements(self):
        """Find out how many output elements there are in total."""

        self.numElem = None   # Placeholder for "not known"
        if self.groupName is None: 
            return
                
        # Break file name into base and running sequence number
        realFileName = os.path.split(fileName)[1]
        fileNameParts = realFileName.split('_')

        if self.groupName[:8] == 'PartType':
            ptIndex = int(self.groupName[8])
        else:
            ptIndex = None   # For checking

        # Deal with particle catalogue files
        if (fileNameParts[0] in ['snap', 'snip', 'partMags',
                                 'eagle_subfind_particles']):
            if verbose: print("Particle catalogue detected!")

            # Need to extract particle index to count (mag --> stars!)
            if self.groupName[:8] == 'PartType':
                self._count_elements_snap(int(self.groupName[8]))
            elif fileNameParts[0] == 'partMags':
                self._count_elements_snap(4)
            else:
                return    # Can't determine element numbers then

        # Deal with subfind catalogue files
        elif len(fileNameParts) >= 2:
            if "_".join(fileNameParts[:3]) == 'eagle_subfind_tab':
                if verbose: print("Subfind catalogue detected!")
                if groupName == 'FOF':
                    self._count_elements_sf_fof()
                elif groupName == 'Subhalo':
                    self._count_elements_sf_subhalo()
                elif groupName == 'IDs':
                    self._count_elements_sf_ids()
                else:
                    return
            
        # In all other cases, don't know how to pre-determine
        # element count, so need to find out from individual files (later).
        else:
            return

        # Final piece: deal with possibility of int32 overflow in numbers
        # (it does happen somewhere...)
        if self.numElem  < 0:
            self.numElem += 4294967296  
        if verbose: print("   ... (detected {:d} elements) ..."
                          .format(self.numElem), flush=True)
    
    def _count_elements_snap(self, pt_index):
        """Count number of particles of specified (int) type."""
        self.numElem = hd.read_attribute(
            self.fileName, 'Header', 'NumPart_Total', require=True)[pt_index]

    def _count_elements_sf_fof(self):
        """Count number of FOF groups in Subfind catalogue."""
        self.numElem = hd.read_attribute(
            self.fileName, 'Header', 'TotNgroups', require=True)

    def _count_elements_sf_subhalo(self):
        """Count number of subhaloes in Subfind catalogue."""
        self.numElem = hd.read_attribute(
            self.fileName, 'Header', 'TotNsubgroups', require=True)

    def _count_elements_sf_ids(self):
        """Count number of particle IDs in Subfind catalogue."""
        self.numElem = hd.read_attribute(
            self.fileName, 'Header', 'TotNids', require=True)

    def _count_files(self, simType='Eagle'):
        """Count number of files to read."""

        self.numFiles = None # Default
        
        if simType == 'Eagle':
            self.numFiles = hd.read_attribute(
                self.fileName, 'Header', 'NumFilesPerSnapshot', require=True)
        elif simType == 'Illustris':
            self.numFiles = hd.read_attribute(
                self.fileName, 'Header', 'NumFiles', require=True)
        else:
            print("Unknown simulation type '{:s}'." .format(simType))
            set_trace()

            
