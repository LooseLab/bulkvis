import h5py
from argparse import ArgumentParser, ArgumentTypeError
import configparser


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
    in_args.add_argument("-i", "--input-dir",
                         help="The directory containing bulk-files for visualisation",
                         type=str,
                         required=True,
                         metavar=""
                         )
    in_args.add_argument("-e", "--export-dir",
                         help="The directory where read-files will be written by bulkvis",
                         type=str,
                         required=True,
                         metavar=""
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
                          help="Path to the config.ini file in your bulkvis installation",
                          default='config.ini',
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
labels = int_labels.union(state_labels)
config = configparser.ConfigParser()
config['data'] = {
    'dir': args.input_dir,
    'out': args.export_dir
}
config['plot_opts'] = {
    'wdg_width': 300,
    'plot_width': 980,
    'plot_height': 800,
    'y_min': 0,
    'y_max': 2200,
    'label_height': 800,
    'upper_cut_off': 2200,
    'lower_cut_off': -1000,
    'output_backend': 'canvas',
}
config['labels'] = {}
for label in labels:
    config['labels'][label] = str(args.flag)
if args.standard:
    for thing in config:
        print('[{t}]'.format(t=thing))
        for t in config[thing]:
            print('{t} = {v}'.format(t=t, v=config[thing][t]))

if args.config:
    with open(args.config, 'w') as configfile:
        config.write(configfile)

print("")
