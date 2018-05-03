from argparse import ArgumentParser, ArgumentTypeError, ArgumentDefaultsHelpFormatter
from collections import defaultdict
import h5py
import pandas as pd
from whale_watch import fuse_reads
import sys
from pathlib import Path
import subprocess
from tqdm import tqdm


def full_path(file):
    return str(Path(file).expanduser().resolve())


def validate_file(file):
    f = Path(full_path(file))
    valid = ['.eps', '.ps', '.tex', '.pdf', '.jpeg', '.tiff', '.png', '.bmp', '.svg', '.wmf']
    if f.suffix not in valid:
        msg = "File doesn't end with one of {v}".format(v=valid)
        raise ArgumentTypeError(msg)
    return str(f)


def trim_series(keep_series, trim_df):
    """
    Given an input of a df, and list of Series, delete Series not in the list
    :param keep_series: list of Series to keep
    :param trim_df: pandas DataFrame
    :return: pandas DataFrame with series removed
    """
    for i, s in trim_df.items():
        if type(s).__name__ == "Series" and s.name not in keep_series:
            trim_df = trim_df.drop(s.name, axis=1)
    return trim_df


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


def prepare_data(seq_sum_df, interval_time, bulkfile, v, sf, run_id):
    seq_sum_df['end_time'] = seq_sum_df['start_time'] + seq_sum_df['duration']
    seq_sum_df = seq_sum_df[seq_sum_df['run_id'] == run_id]
    channels = seq_sum_df['channel'].drop_duplicates(keep='first').values
    dist_dict = defaultdict(list)
    sdist_dict = defaultdict(list)

    for ch in tqdm(channels):
        # format channel
        ch = int(ch)
        ch_str = "Channel_{ch}".format(ch=ch)
        # if v:
        #     print(ch_str)
        # set the paths
        try:
            int_path = bulkfile["IntermediateData"][ch_str]["Reads"]
        except KeyError:
            print("KeyError" + str(ch_str))
        except Exception as e:
            print(type(e).__name__)
            print(e)
            sys.exit()
        try:
            state_path = bulkfile["StateData"][ch_str]["States"]
        except KeyError:
            print("KeyError" + str(ch_str))
        except Exception as e:
            print(type(e).__name__)
            print(e)
            sys.exit()
        # get annotations from int data
        int_fields = ['read_id', 'read_start', 'modal_classification']
        labels, labels_dt = get_annotations(int_path, int_fields, 'modal_classification')
        labels = labels.drop_duplicates(subset="read_id", keep="first")
        labels['read_start'] = labels['read_start'] / sf
        labels['read_id'] = labels['read_id'].str.decode('utf8')
        # get annotations from state data
        state_fields = ['acquisition_raw_index', 'summary_state']
        state_label_df, state_label_dtypes = get_annotations(state_path, state_fields, 'summary_state')
        state_label_df['acquisition_raw_index'] = state_label_df['acquisition_raw_index'] / sf
        state_label_df = state_label_df.rename(
            columns={'acquisition_raw_index': 'read_start', 'summary_state': 'modal_classification'}
        )
        # merge int and state data
        labels = labels.append(state_label_df, ignore_index=True)
        labels.sort_values(by='read_start', ascending=True, inplace=True)
        labels_dt.update(state_label_dtypes)
        labels = labels.fillna(method='ffill')
        # match against channel
        seq_sum_match = seq_sum_df[seq_sum_df['channel'] == ch].dropna()
        for index, row in seq_sum_match.iterrows():
            b = labels[
                labels['read_start'].between(row['end_time'] - interval_time, row['end_time'], inclusive=False)]
            a = labels[
                labels['read_start'].between(row['end_time'], row['end_time'] + interval_time, inclusive=True)]
            for index2, row2 in b.iterrows():
                dist_dict[labels_dt[row2['modal_classification']]].append((row2['read_start'] - row['end_time']))
            for index3, row3 in a.iterrows():
                dist_dict[labels_dt[row3['modal_classification']]].append((row3['read_start'] - row['end_time']))
        # starts
        for index, row in seq_sum_match.iterrows():
            sb = labels[
                labels['read_start'].between(row['start_time'] - interval_time, row['start_time'], inclusive=False)]
            sa = labels[
                labels['read_start'].between(row['start_time'], row['start_time'] + interval_time, inclusive=True)]
            for index2, row2 in sb.iterrows():
                sdist_dict[labels_dt[row2['modal_classification']]].append((row2['read_start'] - row['start_time']))
            for index3, row3 in sa.iterrows():
                sdist_dict[labels_dt[row3['modal_classification']]].append((row3['read_start'] - row['start_time']))

    dist_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dist_dict.items()]))
    sdist_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in sdist_dict.items()]))
    return dist_df, sdist_df


def main():
    args = get_args()
    """
    Use whale_watch to get fused reads data frames
    df2:               just fused reads df
    ss:                sequencing summary df
    chained_reads_ids: read_ids of just fused reads
    """
    df2, ss, chained_read_ids = fuse_reads(args.summary, args.paf, args.distance, 0, False, False)
    keep_cols = ['channel', 'read_id', 'run_id', 'start_time', 'duration']
    # ss2 is un-fused reads
    ss2 = ss[ss['read_id'].isin(chained_read_ids) == False]
    ss2 = ss2.reset_index()
    # ss3 is all to-be-fused reads
    ss3 = ss[ss['read_id'].isin(chained_read_ids) == True]
    # ss4 is 2nd chained read to last chained read
    ss4 = ss3[ss3['read_id'].isin(df2['read_id']) == False]
    #ss5 is all chained reads excluding the last in the chain
    ss5 = ss3[ss3['read_id'].isin(df2['last_read_id']) == False]
    ss5 = ss5.drop_duplicates(subset=keep_cols, keep='first')
    ss6 = ss3[ss3['read_id'].isin(df2['last_read_id'])==True]
    # Open bulkfile
    bulkfile = h5py.File(args.bulk_file, "r")
    sf = int(bulkfile["UniqueGlobalKey"]["context_tags"].attrs["sample_frequency"].decode('utf8'))
    run_id = bulkfile["UniqueGlobalKey"]["tracking_id"].attrs["run_id"].decode('utf8')
    ss4 = trim_series(keep_cols, ss4)
    ss2 = trim_series(keep_cols, ss2)
    ss5 = trim_series(keep_cols, ss5)
    ss6 = trim_series(keep_cols, ss6)
    df2 = trim_series(keep_cols, df2)

    df2['duration'] = pd.to_numeric(df2['duration'], errors='coerce')
    df2['channel'] = pd.to_numeric(df2['channel'], errors='coerce').astype('int64')

    i_dist_df, i_sdist_df = prepare_data(ss2, args.time, bulkfile, args.verbose, sf, run_id)
    k_dist_df, k_sdist_df = prepare_data(ss4, args.time, bulkfile, args.verbose, sf, run_id)
    j_dist_df, j_sdist_df = prepare_data(df2, args.time, bulkfile, args.verbose, sf, run_id)
    l_dist_df, l_sdist_df = prepare_data(ss6, args.time, bulkfile, args.verbose, sf, run_id)
    m_dist_df, m_sdist_df = prepare_data(ss5, args.time, bulkfile, args.verbose, sf, run_id)
    # starts of un-fused reads "Unique read start"
    i_sdist_df.to_csv(args.A, sep=",", header=True)
    # unique read ends
    i_dist_df.to_csv(args.B, sep=",", header=True)
    # start of just fused reads
    j_sdist_df.to_csv(args.C, sep=",", header=True)
    # starts of reads to be fused "Internal read starts"
    k_sdist_df.to_csv(args.D, sep=",", header=True)
    # internal read ends
    m_dist_df.to_csv(args.E, sep=",", header=True)
    # last fused read end
    l_dist_df.to_csv(args.F, sep=",", header=True)

    if not args.no_generate_plot:
        # subprocess to R, generate plot
        R_script = str(Path(Path(__file__).resolve().parent / 'whale.R'))
        cmd = "Rscript {r} {a} {b} {c} {d} {e} {f} {o} {run}".format(r=R_script,
                                                                     a=args.A,
                                                                     b=args.B,
                                                                     c=args.C,
                                                                     d=args.D,
                                                                     e=args.E,
                                                                     f=args.F,
                                                                     o=args.out,
                                                                     run=run_id
                                                                     )
        proc = subprocess.Popen(cmd, shell=True)
        proc.communicate()
        if proc.returncode == 0:
            print("plot saved as {f}".format(f=Path(args.out).name))
        else:
            print(proc.returncode)
    return


def get_args():
    parser = ArgumentParser(
        description="""Parse sequencing_summary.txt, .paf, and 
                    bulk fast5 files to generate CSV files 
                    containing the distributions of MinKNOW 
                    events around read starts and ends. These 
                    are divided into unique reads, split reads 
                    and internal reads. The R script, whale.R, 
                    is called to generate the plot; this requires
                    the packages: ggplot2, tidyr, and dplyr.
                    Note: of the MinKNOW classifications only above, 
                    adapter, pore, transition, unblocking, and 
                    unclassified are included.""",
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False)
    general = parser.add_argument_group(
        title='General options')
    general.add_argument("-h", "--help",
                         action="help",
                         help="Show this help and exit"
                         )
    general.add_argument("-d", "--distance",
                         help='''Specify the maximum distance, in bases, between consecutive mappings.
                              This is the difference between \'Target Start\' and \'Target End\' in 
                              a paf file''',
                         type=int,
                         default=10000,
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
                         help="An ONT bulk fast5 file containing raw signal",
                         type=str,
                         required=True,
                         )
    in_args.add_argument("-s", "--summary",
                         help="A sequencing summary file generated by albacore",
                         type=str,
                         required=True,
                         )
    in_args.add_argument("-p", "--paf",
                         help="A paf file generated by minimap2",
                         type=str,
                         required=True,
                         )
    in_args.add_argument("-t", "--time",
                         help='''+/- time around a strand event in seconds''',
                         type=int,
                         required=True,
                         default=10,
                         )
    out_args = parser.add_argument_group(
        title='Output files'
    )
    out_args.add_argument("--no-generate-plot",
                          help="If set, do not generate density plot",
                          action="store_true",
                          default=False,
                          )
    out_args.add_argument("-A",
                         help="CSV of MinKNOW events occurring before and after correctly called read starts",
                         type=full_path,
                         default='unique_read_start.csv',
                         required=False,
                         )
    out_args.add_argument("-B",
                         help="CSV of MinKNOW events occurring before and after correctly called read ends",
                         type=full_path,
                         default='unique_read_end.csv',
                         required=False,
                         )
    out_args.add_argument("-C",
                         help="""CSV of MinKNOW events occurring before and after 
                         the start of the first incorrectly split read in a group""",
                         type=full_path,
                         default='split_read_start.csv',
                         required=False,
                         )
    out_args.add_argument("-D",
                         help="""CSV of MinKNOW events occurring before and after incorrectly 
                         called read starts, within a group of incorrectly split reads""",
                         type=full_path,
                         default='internal_read_start.csv',
                         required=False,
                         )
    out_args.add_argument("-E",
                          help="""CSV of MinKNOW events occurring before and after incorrectly 
                          called read ends, within a group of incorrectly split reads""",
                          type=full_path,
                          default="internal_read_end.csv",
                          required=False,
                           )
    out_args.add_argument("-F",
                         help="""CSV of MinKNOW events occurring before and after the end 
                         of the first incorrectly split read in a group""",
                         type=full_path,
                         default='split_read_end.csv',
                         required=False,
                         )
    out_args.add_argument("--out",
                          help='''Specify the output filename for the plot.
                          File extension must be one of [.eps, .ps, .tex, .pdf, .jpeg, .tiff, .png, .bmp, .svg, .wmf]''',
                          type=validate_file,
                          default='classification_count.pdf',
                          required=False,
                           )
    return parser.parse_args()


if __name__ == '__main__':
    main()
