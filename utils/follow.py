import h5py
import pandas as pd
from argparse import ArgumentParser

def get_annotations(path, fields, enum_field):
    data_labels = {}
    for field in fields:
        data_labels[field] = path[field]
    data_dtypes = {}
    if h5py.check_dtype(enum=path.dtype[enum_field]):
        dataset_dtype = h5py.check_dtype(enum=path.dtype[enum_field])
        # data_dtype may lose some dataset dtypes there are duplicates of 'v'
        data_dtypes = {v: k for k, v in dataset_dtype.items()}
    labels_df = pd.DataFrame(data=data_labels)
    return labels_df, data_dtypes, dataset_dtype


def main():
    args = get_args()
    bulkfile = h5py.File(args.bulk_file, "r")
    sf = int(bulkfile["UniqueGlobalKey"]["context_tags"].attrs["sample_frequency"].decode('utf8'))
    run_id = bulkfile["UniqueGlobalKey"]["tracking_id"].attrs["run_id"].decode('utf8')

    ss_fields = ['channel', 'read_id', 'run_id']
    seq_sum_df = pd.read_csv(args.summary, sep='\t', usecols=ss_fields)
    seq_sum_df = seq_sum_df[seq_sum_df['run_id'] == run_id]
    seq_sum_df['ch_str'] = "Channel_" + seq_sum_df['channel'].map(str)

    q_keys = pd.DataFrame(bulkfile['IntermediateData']['Channel_1']['Reads'][()])
    q_keys = q_keys.keys()

    labels, label_dtypes, rev_dtypes = get_annotations(
        bulkfile['IntermediateData']['Channel_1']['Reads'],
        '',
        'modal_classification'
    )
    chans = seq_sum_df['ch_str'].drop_duplicates(keep='first').values
    result = pd.DataFrame(columns=q_keys)
    result_arr = []
    for ch in chans:
        try:
            df = pd.DataFrame(bulkfile['IntermediateData'][ch]['Reads'][()]).drop_duplicates(
                subset="read_id",
                keep="first"
            )
            result_arr.append(df)
        except KeyError:
            pass
        except Exception as e:
            print(type(e).__name__)

    result = pd.concat(result_arr)
    
    fields = [('P1', 1), ('P2', 2), ('P3', 3), ('P4', 4), ('P5', 5),
              ('F1', -1), ('F2', -2), ('F3', -3), ('F4', -4), ('F5', -5)]
    
    result['read_id'] = result['read_id'].str.decode('utf8')

    for f, s in fields:
        result[f] = result['modal_classification'].shift(s)

    dicta = {}
    dictb = {}
    result = result[result['read_id'].isin(seq_sum_df['read_id'].values)]
    for f, s in fields:
        dicta[f] = result[f].value_counts().to_dict()
        dictb['{f}_1'.format(f=f)] = {}
        for k, v in dicta[f].items():
            dictb['{f}_1'.format(f=f)][label_dtypes[int(k)]] = v

    print(dictb)

def get_args():
    parser = ArgumentParser(
        description="""This does something...""",
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
    in_args.add_argument("-b", "--bulk-file",
                         help="An ONT bulk-file containing raw signal",
                         type=str,
                         default='',
                         required=True,
                         metavar=''
                         )
    in_args.add_argument("-s", "--summary",
                         help="An ONT sequencing summary file",
                         type=str,
                         default='',
                         required=True,
                         metavar=''
                         )
    return parser.parse_args()


if __name__ == '__main__':
    main()