# Blank file

import os
from pdb import set_trace

class SplitFile:
    """Class to read data sets that are split over multiple files.
    """
        
    def __init__(self, fileName, groupName, ptype=None):
                
        # First: check if the file exists...!
        if not os.path.isfile(fileName):
            print("-----------------------------------------------------")
            print("Error [" + fileName + "]")
            print("This file does not exist. Please supply correct name!")
            print("-----------------------------------------------------")
            return 
        self.fileName = fileName
        
        # Set up numerical particle type, if appropriate
        if ptype is not None:
            self._decode_ptype(ptype)
        else:
            self.groupName = groupName
            if groupName[:8] == 'PartType':
                self.groupIndex = int(groupName[8])
        
        if not silent:
            print("Preparing reading from '{:s}'..."
                  .format(groupName.upper()))

        # Find the total number of elements in the output array
        self._count_elements()


    def _count_elements(self):
        """Find out how many output elements there are in total."""
        
        # Break file name into base and running sequence number
        realFileName = os.path.split(fileName)[1]
        fileNameParts = realFileName.split('_')
        
        self.typeCode = 0 # Set as 'unknown' by default
    
        # ------- Determine type of file to be read -----------

        if (fileNameParts[0] == 'snip'):
            if not silent: print("SNIPSHOT detected!")
            self._count_elements_snip()
        elif (fileNameParts[0] == 'snap'):
            if not silent: print("SNAPSHOT detected!")
            self._count_elements_snap()
        elif (type_parts[0] == 'partMags'):
            if not silent: print("PARTMAGS detected!")
            type_code = 2
            self.groupIndex = 4   # Need to set here.
        astro = False   # None of these data need astro corrections...
    elif(type_parts[0] in ['rockstar', 'Rockstar']):
        if not silent: print("ROCKSTAR detected!")
        type_code = 5
    elif len(type_parts) >= 2:
        if ("_".join(type_parts[0:3]) == 'eagle_subfind_tab'):
            if not silent: print("EAGLE_SUBFIND_TAB detected!")
            type_code = 3
        elif ("_".join(type_parts[0:3]) == 'eagle_subfind_particles'):
            if not silent: print("EAGLE_SUBFIND_PARTICLES detected!")
            type_code = 4



            
    def _decode_ptype(self, ptype):
        """Identify supplied particle type"""
        
        if isinstance(ptype, int):
            self.groupName = 'PartType{:d}' .format(ptype)
            self.groupIndex = ptype
        elif isinstance(ptype, str):
            if ptype.upper() == 'GAS':
                self.groupName = 'PartType0'
                self.groupIndex = 0
            elif ptype.upper() == 'DM':
                self.groupName = 'PartType1'
                self.groupIndex = 1
            elif ptype.upper() in ['STARS', 'STAR']:
                self.groupName = 'PartType4'
                self.groupIndex = 4
            elif ptype.upper() in ['BH', 'BHS', 'BLACKHOLE', 'BLACKHOLES',
                                   'BLACK_HOLES', 'BLACK_HOLE']:
                self.groupName = 'PartType5'
                self.groupIndex = 5
            else:
                print("Unrecognized particle type name '{:s}'"
                      .format(ptype))
                set_trace()
        else:
            print("Unrecognized particle type format")
            set_trace()
