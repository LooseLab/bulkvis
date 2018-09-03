import h5py
from argparse import ArgumentParser, ArgumentTypeError
from pathlib import Path
import configparser


def full_path(file):
    q = Path(file).expanduser().resolve()
    if q.exists() and q.is_dir():
        return str(q)
    else:
        msg = "Path '{v}' doesn't exist or isn't a directory".format(v=file)
        raise ArgumentTypeError(msg)


def get_args():
    parser = ArgumentParser(
        description="""Generate a configuration file required for bulkvis to run""",
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
                         required=True,
                         metavar=''
                         )
    in_args.add_argument("-i", "--input-dir",
                         help="The path to tbe folder containing bulk-files for visualisation",
                         type=full_path,
                         required=True,
                         metavar=""
                         )
    in_args.add_argument("-e", "--export-dir",
                         help="The path to tbe folder where read-files will be written by bulkvis",
                         type=full_path,
                         required=True,
                         metavar=""
                         )
    in_args.add_argument("-m", "--map-dir",
                         help="The path to tbe folder where map files, created by gen_bmf.py, are stored",
                         type=full_path,
                         required=True,
                         metavar=""
                         )
    out_args = parser.add_argument_group(
        title='Output'
    )
    out_args.add_argument('-c', '--config',
                          help="Path to the config.ini file in your bulkvis installation",
                          default='config.ini',
                          type=str,
                          metavar=''
                          )
    return parser.parse_args()


def get_annotations(path, enum_field):
    data_dtypes = {}
    if h5py.check_dtype(enum=path.dtype[enum_field]):
        data_dtypes = h5py.check_dtype(enum=path.dtype[enum_field])
    return data_dtypes

def main():
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
        'out': args.export_dir,
        'map': args.map_dir
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
    default_opts = {'adapter', 'mux_uncertain', 'user1', 'pore', 'strand', 'transition', 'unblocking', 'above'}
    for label in labels:
        if label in default_opts:
            config['labels'][label] = 'True'
        else:
            config['labels'][label] = 'False'

    if args.config:
        with open(args.config, 'w') as configfile:
            config.write(configfile)

if __name__ == '__main__':
    main()
