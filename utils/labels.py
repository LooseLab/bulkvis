import h5py
from argparse import ArgumentParser, ArgumentTypeError
import configparser
import sys


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
                                  overwrite any current label configuration in this file.''',
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

def get_annotations(path, enum_field):
    data_dtypes = {}
    if h5py.check_dtype(enum=path.dtype[enum_field]):
        data_dtypes = h5py.check_dtype(enum=path.dtype[enum_field])
    return data_dtypes


args = get_args()
bulkfile = h5py.File(args.bulkfile, 'r')
int_data_path = bulkfile["IntermediateData"]
state_data_path = bulkfile["StateData"]

for i, channel in enumerate(int_data_path):
    if i == 0:
        path = int_data_path[channel]["Reads"]
        int_labels = get_annotations(path, 'modal_classification')
    break

for i, channel in enumerate(state_data_path):
    if i == 0:
        path = state_data_path[channel]["States"]
        state_labels = get_annotations(path, 'summary_state')
    break
int_labels = set(int_labels.keys())
state_labels = set(state_labels.keys())
# concat sets here
# Move STDOUT to end of file
if args.standard:
    for label in int_labels:
        print(label)
    for label in state_labels:
        print(label)

if args.config:
    """
    Open config file
    check and preserve all sections other than 'labels'
    clear labels
    write from sets to labels using -f for state
    """