"""
DOCSTRING: set_config.py
This script generates and populates a config.ini file
it requires a bulkfile as input material to set available
labels.
"""

from argparse import ArgumentParser, ArgumentTypeError
import h5py
import configparser


"""
general action:
is there a config file
"""



def get_args():
    parser = ArgumentParser(
        description="""Get the label options for a bulkfile for copying into a config file""",
        add_help=False)
    general = parser.add_argument_group(
        title='General options')
    general.add_argument("-h", "--help",
                         action="help",
                         help="Show this help and exit"
                         )
    in_args = parser.add_argument_group(
        title='Input sources'
    )
    in_args.add_argument("-b", "--bulkfile",
                         help="A bulk-fast5 file to get labels from",
                         type=str,
                         default="",
                         required=True,
                         metavar=''
                         )
    out_args = parser.add_argument_group(
        title='Output'
    )
    out_args.add_argument('-s', '--standard',
                          help='''Output to STDOUT''',
                          action="store_true",
                          default=False,
                          )
    out_args.add_argument('-c', '--config',
                          help='''Specify a config file to write the labels to. This will 
                                  overwrite any current label configuration in this file.
                                  If no config file is specified one will be saved in the 
                                  current terminal directory.''',
                          type=str,
                          metavar=''
                          )
    out_args.add_argument('-f', '--flag',
                          help='''Set label flag as, defaults to True''',
                          type=str2bool,
                          default=True,
                          metavar=''
                          )
    return parser.parse_args()

def str2bool(v):
    if v.lower() in ('true', 't', 'yes', 'y', '1'):
        return True
    elif v.lower() in ('false', 'f', 'no', 'n', '0'):
        return False
    else:
        return ArgumentTypeError('Boolean value expected.')