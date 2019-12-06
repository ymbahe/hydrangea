# Blank file

def clone_dir(dir, loc = 'freya'):
    """
    Construct a 'clone' of a specified directory structure.

    This is used to enable storing experimental/development/specialized outputs
    in separate locations from the 'main' simulation repository, but with an
    exactly analogous internal structure.

    Args:
        dir (string): The full path of the directory to clone.
        loc (string, optional): The branch to clone the directory to. Can be
            one of 'freya' or 'virgo'.
    
    Returns:
        string: modified directory name in the specified branch.
    """
        
    dir_parts = dir.split('/')

    num_Hydrangea = dir_parts.count('Hydrangea')
    if dir_parts.count('Hydrangea') == 0:
        print("The input path '" + dir + "' does not seem to contain a directory called 'Hydrangea'...")
        sys.exit(44)

    last_of_pre = dir_parts.index('Hydrangea')
    first_special = last_of_pre + 1

    special_part = '/'.join(dir_parts[first_special:])

    if loc == 'freya':
        return "/freya/ptmp/mpa/ybahe/HYDRANGEA/ANALYSIS/" + special_part
    elif loc == 'virgo':
        return "/virgo/scratch/ybahe/HYDRANGEA/ANALYSIS/" + special_part
    else:
        print("Do not understand requested redirect site '" + loc + "'")
        set_trace()

def get_snepshot_indices(rundir, list='basic'):
    """
    Extract type, number, and aexp for snepshots from a specified list.
    """
    
    snepdir = rundir + '/sneplists/'
    fileName = snepdir + list + '.dat'

    data = ascii.read(fileName)
    
    rootIndex = np.array(data['rootIndex'])
    aexp = np.array(data['aexp'])
    sourceType = np.array(data['sourceType'])
    sourceNum = np.array(data['sourceNum'])
    
    return rootIndex, aexp, sourceType, sourceNum
    

def snap_times(conv = None, list = None):
    """
    Return the times of all Hydrangea snapshots.

    By default, the expansion factors of the 30 snapshots are returned.
    Optionally, two arguments can be provided:

    Args:
        conv (string, 'zred' or 'age'): convert the expansion factors to 
            redshift or age of the Universe (as appropriate).
        list (string): name (without directories) of the output file to read.
    """
    
    if list is None:
        snaptimes_file = currDir + '/hydrangea/OutputLists/hydrangea_snapshots_plus.dat'
    else:
        snaptimes_file = currDir + '/hydrangea/OutputLists/' + list + '.dat'

    snap_times = np.array(ascii.read(snaptimes_file, format = 'no_header')['col1'])

    if conv is None:
        return snap_times

    elif conv == 'zred':
        return 1/snap_times - 1

    elif conv == 'age':
        return Planck13.age(1/snap_times - 1).to(u.Gyr).value

    else:
        print("I do not know what you mean by '" + conv + "'...")
        set_trace()

    
def test():
    print("Can you see me still?")
