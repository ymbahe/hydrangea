# Blank file

from astropy.io import ascii
from astropy.cosmology import Planck13
import numpy as np
import os
from pdb import set_trace()

def get_snepshot_indices(rundir, list='basic'):
    """
    Extract type, number, and aexp for snepshots from a specified list.

    Parameters
    ----------
    rundir : str
        The base directory of the simulation.
    list : str, optional
        The snepshot list to extract (default: 'basic')

    Returns
    -------
    rootIndex : np.array(int)
        The root indices of all snepshots in specified list.
    aExp : np.array(float):
        The expansion factor of all snepshots in specified list.
    sourceType : np.array(string)
        For each snepshot, whether it is a 'snap' or 'snip'.
    sourceNum : np.array(int)
        Indices of snepshot in its category (i.e. snap/snipshot).
    """
    
    snepdir = rundir + '/sneplists/'
    fileName = snepdir + list + '.dat'

    data = ascii.read(fileName)
    
    rootIndex = np.array(data['rootIndex'])
    aExp = np.array(data['aexp'])
    sourceType = np.array(data['sourceType'])
    sourceNum = np.array(data['sourceNum'])
    
    return rootIndex, aexp, sourceType, sourceNum
    
def snep_times(timeType='aexp', list='allsnaps'):
    """Return the times of all snepshots in a particular list.

    By default, the expansion factors of the 30 snapshots are returned.

    Parameters
    ----------
    timeType : string
        Specified the 'time' flavour to be returned. Options are:
            'aexp' [default]: Expansion factor
            'zred': Redshift
            'age': Age of the Universe [Gyr]
            'lbt': Lookback time from z = 0 [Gyr]
    list : string
        The snepshot list for which to load times. Options are:
            'z0_only': 1 snapshot at z = 0
            'regsnaps': 28 regular snapshots
            'allsnaps' [default]: 30 snapshots
            'basic': 109 snepshots at Delta_t = 125 Myr
            'default_long': basic + filler snipshots to Delta_t = 25 Myr
            'short_movie': 217 snepshots at Delta_t = 62.5 Myr
            'full_movie': 1081 snepshots at Delta_t = 12.5 Myr

    Returns
    -------
    np.array(float)
        The desired snepshot times.

    Note
    ----
    This function retrieves the target time of the snepshot, which may 
    deviate slightly from the the actual output. Use read_snepshot_time() 
    for the latter.
    
    """

    currDir = os.path.dirname(os.path.realpath(__file__)) + "/"
    snaptimes_file = currDir + 'OutputLists/' + list + '.dat'

    aExp = np.array(ascii.read(snaptimes_file, format = 'no_header')['col1'])

    if timeType == 'aexp':
        return aExp
    elif timeType == 'zred':
        return 1/aExp - 1
    elif timeType == 'age':
        return Planck13.age(1/aExp - 1).to(u.Gyr).value
    elif timeType == 'lbt':
        return Planck13.lookback_time(1/aExp - 1).to(u.Gyr).value
    else:
        print("I do not know what you mean by '" + timeType + "'...")
        set_trace()
