import sys

__all__=["pdb","pdb_viewer","pdb_basic","molecule","chain","residue","atom"]

__version__ = "2.0.0"

def Minium_version_required(ver_str):
    if [int(x) for x in __version__.split(".")] < [int(x) for x in ver_str.split(".")]:
        sys.exit("Error: Minimum requirement of molar module version is %s.\
    The installed version is %s\
    Please update your molar module from: https://github.com/alinar/Molar" % (ver_str , molar.__version__) )
    else:
        return True