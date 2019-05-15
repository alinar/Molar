__all__=["pdb","pdb_viewer","pdb_basic","molecule","chain","residue","atom"]
__version__ = "2.1.4"

import sys

def Minium_version_required(ver_str):
    if [int(x) for x in __version__.split(".")] < [int(x) for x in ver_str.split(".")]:
        sys.exit("Error: The minimum required version of \"molar\" module is %s.\n\
The installed version is %s\n\
Please update your molar module from: https://github.com/alinar/Molar\n\
Or use the command:\n\
sudo pip install git+https://github.com/alinar/Molar\n\
\nThe path of the module in use is: %s" % ( ver_str,__version__,__path__[0] ) )
    else:
        return True
