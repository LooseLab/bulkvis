import h5py
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from argparse import ArgumentParser
import sys


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
    return labels_df, data_dtypes


def prepare_data(summary_file, interval_time, bulkfile, v, sf, run_id):
    ss_fields = ['channel', 'read_id', 'run_id', 'start_time']
    seq_sum_df = pd.read_csv(summary_file, sep='\t', usecols=ss_fields)
    seq_sum_df = seq_sum_df[seq_sum_df['run_id'] == run_id]
    channels = seq_sum_df['channel'].drop_duplicates(keep='first').values
    # interval_time = args.time
    b_dict = defaultdict(lambda: defaultdict(int))
    a_dict = defaultdict(lambda: defaultdict(int))
    qi = 0
    for ch in channels:
        qi += 1
        ch_str = "Channel_{ch}".format(ch=ch)
        if v:
            print(ch_str)
        try:
            int_path = bulkfile["IntermediateData"][ch_str]["Reads"]
            state_path = bulkfile["StateData"][ch_str]["States"]
        except KeyError:
            print("KeyError" + str(ch_str))
        except Exception as e:
            print(type(e).__name__)
            print(e)
            sys.exit()
        int_fields = ['read_id', 'read_start', 'modal_classification']
        labels, labels_dt = get_annotations(int_path, int_fields, 'modal_classification')
        labels = labels.drop_duplicates(subset="read_id", keep="first")
        labels['read_start'] = labels['read_start'] / sf
        labels['read_id'] = labels['read_id'].str.decode('utf8')

        state_fields = ['acquisition_raw_index', 'summary_state']
        state_label_df, state_label_dtypes = get_annotations(state_path, state_fields, 'summary_state')
        state_label_df['acquisition_raw_index'] = state_label_df['acquisition_raw_index'] / sf
        state_label_df = state_label_df.rename(
            columns={'acquisition_raw_index': 'read_start', 'summary_state': 'modal_classification'}
        )
        labels = labels.append(state_label_df, ignore_index=True)
        labels.sort_values(by='read_start', ascending=True, inplace=True)
        labels_dt.update(state_label_dtypes)
        labels = labels.fillna(method='ffill')
        seq_sum_match = seq_sum_df[seq_sum_df['channel'] == ch].dropna()
        for index, row in seq_sum_match.iterrows():
            b = labels[
                labels['read_start'].between(row['start_time'] - interval_time, row['start_time'], inclusive=False)]
            a = labels[
                labels['read_start'].between(row['start_time'], row['start_time'] + interval_time, inclusive=False)]
            b_count = 0
            for index2, row2 in b.iterrows():
                b_count -= 1
                if b_count >= -5:
                    b_dict[labels_dt[row2['modal_classification']]][b_count] += 1
            a_count = 0
            for index3, row3 in a.iterrows():
                a_count += 1
                if a_count <= 5:
                    a_dict[labels_dt[row3['modal_classification']]][a_count] += 1

    merged_dict = defaultdict(dict)

    merged_dict.update(a_dict)
    for key, nested_dict in b_dict.items():
        merged_dict[key].update(nested_dict)

    df = pd.DataFrame(merged_dict)
    df = df.fillna(0)
    df['rowsum'] = df.sum(axis=1)
    df = df.div(df['rowsum'], axis=0)
    df.drop('rowsum', axis=1, inplace=True)
    if v:
        print(df)
    return df


def main():
    args = get_args()
    bulkfile = h5py.File(args.bulk_file, "r")
    sf = int(bulkfile["UniqueGlobalKey"]["context_tags"].attrs["sample_frequency"].decode('utf8'))
    run_id = bulkfile["UniqueGlobalKey"]["tracking_id"].attrs["run_id"].decode('utf8')

    i_df = prepare_data(args.summary_1, args.time, bulkfile, args.verbose, sf, run_id)
    i_df.to_csv("just_fused_prop.csv", sep=",", header=True)
    j_df = prepare_data(args.summary_2, args.time, bulkfile, args.verbose, sf, run_id)
    j_df.to_csv("to_be_fused_prop.csv", sep=",", header=True)
    k_df = prepare_data(args.summary_3, args.time, bulkfile, args.verbose, sf, run_id)
    k_df.to_csv("un_fused_prop.csv", sep=",", header=True)


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
    general.add_argument("-V", "--verbose",
                         help='''Print verbose output to terminal''',
                         action="store_true",
                         default=False,
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
    in_args.add_argument("-s1", "--summary-1",
                         help="An ONT sequencing summary file",
                         type=str,
                         default='',
                         required=True,
                         metavar=''
                         )
    in_args.add_argument("-s2", "--summary-2",
                         help="An ONT sequencing summary file",
                         type=str,
                         default='',
                         required=True,
                         metavar=''
                         )
    in_args.add_argument("-s3", "--summary-3",
                         help="An ONT sequencing summary file",
                         type=str,
                         default='',
                         required=True,
                         metavar=''
                         )
    general.add_argument("-t", "--time",
                         help='''plus/minus time around a strand event''',
                         type=int,
                         required=True,
                         default=5,
                         metavar=''
                         )
    return parser.parse_args()


if __name__ == '__main__':
    main()