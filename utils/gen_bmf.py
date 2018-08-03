import pandas as pd
import numpy as np
from argparse import ArgumentParser
from pathlib import Path
import warnings


def main():
    warnings.warn('This is temporary and will be removed without notice')
    args = get_args()
    paf_path = args.paf
    col_names = ['Qname', 'Strand', 'Tname', 'Tstart', 'Tend', 'alignment_type']
    pf = pd.read_csv(paf_path, sep='\t', header=None, names=col_names, usecols=[0, 4, 5, 7, 8, 12])
    pf = pf[pf['alignment_type'] == 'tp:A:P']
    pf = pf.drop_duplicates()

    cols = ['read_id', 'run_id', 'channel', 'start_time', 'duration']
    ss = pd.read_csv(args.summary, sep='\t', usecols=cols)

    df = pd.merge(ss, pf, left_on='read_id', right_on='Qname', how='outer')
    df = df.dropna()

    df['end_time'] = df['start_time'] + df['duration']
    df['start_mapping'] = np.where(df['Strand'] == '+',
                                   df[['Tstart', 'Tend']].min(axis=1),
                                   df[['Tstart', 'Tend']].max(axis=1)
                                   )
    df['end_mapping'] = np.where(df['Strand'] == '-',
                                 df[['Tstart', 'Tend']].min(axis=1),
                                 df[['Tstart', 'Tend']].max(axis=1)
                                 )
    df['sm'] = df['start_mapping'].astype('int64').map('{0:,d}'.format)
    df['em'] = df['end_mapping'].astype('int64').map('{0:,d}'.format)
    df['label'] = (df['Tname'].astype('str') + ": " +
                   df['sm'].astype('str') + " - " +
                   df['em'].astype('str')
                   )

    df = df.rename(columns={'Tname': 'target_name', 'Strand': 'strand'})
    # export csv as some filetype...
    header = ['run_id', 'read_id', 'channel', 'start_time',
              'end_time', 'target_name', 'strand',
              'start_mapping', 'end_mapping', 'label']
    df.to_csv(args.bmf, sep="\t", header=True, columns=header, index=False)


def full_path(file):
    return str(Path(file).expanduser().resolve())


def get_args():
    parser = ArgumentParser(
        description="""Parse sequencing_summary.txt files 
                       and .paf files to format mapping info for bulkvis""",
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
    in_args.add_argument("-s", "--summary",
                         help="A sequencing summary file generated by albacore",
                         type=full_path,
                         default="",
                         required=True,
                         metavar=''
                         )
    in_args.add_argument("-p", "--paf",
                         help="A paf file generated by minimap2",
                         type=full_path,
                         default='',
                         required=True,
                         metavar=''
                         )
    out_args = parser.add_argument_group(
        title='Output files'
    )
    out_args.add_argument('--bmf',
                          help='''Specify the output file name, default is \'mapping.bmf\'''',
                          type=full_path,
                          default='mappings.bmf',
                          metavar=''
                          )
    return parser.parse_args()


if __name__ == '__main__':
    main()

