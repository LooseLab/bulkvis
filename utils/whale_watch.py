from argparse import ArgumentParser
from collections import OrderedDict
import numpy as np
import pandas as pd

def main():
    args = get_args()
    if args.debug:
        debug(args)

    sequencing_summary = args.summary
    ss_fields = ['channel', 'start_time', 'duration', 'run_id', 'read_id', 'sequence_length_template', 'filename']
    ss = pd.read_csv(sequencing_summary, sep='\t', usecols=ss_fields)
    ss = ss.sort_values(by=['channel', 'run_id', 'start_time'])
    ss['diff'] = ss['start_time'].shift(-1) - (ss['start_time'] + ss['duration'])
    # quickly put next read_id on as another column
    ss['next_read_id'] = ss['read_id'].shift(-1)
    ss['next_start_time'] = ss['start_time'].shift(-1)
    ss['end_time'] = ss['start_time'] + ss['duration']
    ss['next_end'] = ss['next_start_time'] + ss['duration'].shift(-1)
    ss['next_sequence_length_template'] = ss['sequence_length_template'].shift(-1)
    ss['combined_length'] = ss['sequence_length_template'] + ss['next_sequence_length_template']
    ss['next_sequence_length_template'] = ss['next_sequence_length_template'].fillna(0).astype('int64')
    ss['combined_length'] = ss['combined_length'].fillna(0).astype('int64')

    stats = OrderedDict()
    stats_n50 = OrderedDict()
    stats['Original read count:'] = len(ss)
    stats_n50['Original N50:'] = n50(ss['sequence_length_template'])

    paf_file = args.paf
    pf = pd.read_csv(paf_file, sep='\t', header=None, usecols=[0, 4, 5, 7, 8])
    pf.columns = ['Qname', 'Strand', 'Tname', 'Tstart', 'Tend']
    """
    0		Query sequence name
    1x		Query sequence length
    2x		Query start coordinate (0-based)
    3x		Query end coordinate (0-based)
    4		`+' if query and target on the same strand; `-' if opposite
    5		Target sequence name (~chr)
    6x		Target sequence length
    7		Target start coordinate on the original strand
    8		Target end coordinate on the original strand
    9x		Number of matching bases in the mapping
    10x		Number bases, including gaps, in the mapping
    11x		Mapping quality (0-255 with 255 for missing)
    """
    df = pd.merge(ss, pf, left_on='read_id', right_on='Qname', how='outer')
    df2 = pd.merge(df, pf, left_on='next_read_id', right_on='Qname', how='outer', suffixes=("_A", "_B"))
    df2 = df2.dropna().reset_index()

    not_qname = df2['Qname_A'] != df2['Qname_B']
    yes_strand = df2['Strand_A'] == df2['Strand_B']
    yes_tname = df2['Tname_A'] == df2['Tname_B']

    df2 = df2[not_qname & yes_strand & yes_tname]
    df2['match_distance'] = np.where(
        df2['Strand_A'] == '+',             # condition
        df2['Tstart_B'] - df2['Tend_A'],    # True
        df2['Tstart_A'] - df2['Tend_B']     # False
    )
    df2 = df2[df2['match_distance'] > 0]
    df2 = df2[df2['match_distance'] < args.distance]
    df2 = df2.drop_duplicates(subset=['channel', 'start_time', 'duration',
                                      'next_start_time', 'diff', 'read_id',
                                      'next_read_id', 'sequence_length_template',
                                      'next_sequence_length_template', 'combined_length'],
                              keep='first'
                              )
    # separate df into read groups and set index to cs to allow grouping
    cond_1 = df2['next_read_id'] == df2['read_id'].shift(-1)
    cond_2 = df2['read_id'] == df2['next_read_id'].shift(-1)
    df2['COND'] = np.where(cond_1 | cond_2, True, False)
    df2['W'] = np.where(df2['COND'].shift(1) == False, 1, 0)
    df2['cs'] = df2['W'].cumsum()
    df2 = df2.set_index('cs')
    # group and concatenate read ids
    df2['all_but_last'] = df2.groupby(level='cs')['read_id'].apply('|'.join)
    df2['last_read_id'] = df2.groupby(level='cs')['next_read_id'].last()
    df2['cat_read_id'] = df2['all_but_last'] + "|" + df2['last_read_id']
    # group and combine length
    df2['combined_length'] = df2.groupby(level='cs')['sequence_length_template'].sum()
    df2['last_length'] = df2.groupby(level='cs')['next_sequence_length_template'].last()
    df2['combined_length'] = df2['combined_length'] + df2['last_length']
    # group and add start of mapping and end of mapping
    # take max/min for end/start match from grouped value list
    df2['int_start_match'] = df2.groupby(level='cs')['Tstart_A'].first()
    df2['int_end_match'] = df2.groupby(level='cs')['Tend_B'].last()
    df2['start_match'] = np.minimum(df2['int_start_match'], df2['int_end_match'])
    df2['end_match'] = np.maximum(df2['int_start_match'], df2['int_end_match'])

    # group and add start and end times
    df2['start_time'] = df2.groupby(level='cs')['start_time'].first()
    df2['next_end'] = df2.groupby(level='cs')['next_end'].last()
    # add the duration (time between start and end)
    df2['duration'] = df2['next_end'] - df2['start_time']
    # format and add coordinates
    df2['stime_floor'] = np.floor(df2['start_time']).astype('int64').astype('str')
    df2['etime_ceil'] = np.ceil(df2['next_end']).astype('int64').astype('str')
    df2['channel'] = df2['channel'].astype('int64').astype('str')
    df2['combined_length'] = df2['combined_length'].astype('int64')
    df2['start_match'] = df2['start_match'].astype('int64').astype('str')
    df2['end_match'] = df2['end_match'].astype('int64').astype('str')
    df2['duration'] = df2['duration'].map('{:,.5f}'.format)
    df2['coords'] = df2['channel'] + ":" + df2['stime_floor'] + "-" + df2['etime_ceil']
    # split read-file names to give the bulk-file names
    df2['filename'] = df2['filename'].str.split('_read', 1).str[0]
    df2.rename(columns={'Tname_A': 'target_name', 'Strand_A': 'strand'}, inplace=True)
    # fused_read_ids is a pd.Series of all fused reads
    fused_read_ids = pd.concat([df2['read_id'], df2['next_read_id']])

    df2['count'] = df2.groupby(level='cs')['read_id'].size() + 1

    # remove duplicate entries from df2
    df2 = df2.drop_duplicates(subset=['coords', 'channel', 'start_time',
                                      'duration', 'combined_length',
                                      'start_match', 'end_match', 'cat_read_id'],
                              keep='first')
    fused_read_ids = fused_read_ids.unique()

    # un_fused_df contains reads that are correctly split
    un_fused_df = ss[ss['read_id'].isin(fused_read_ids) == False]
    un_fused_df = un_fused_df.reset_index()
    # split_df is reads that have false starts (i.e 2->N)
    split_df = ss[ss['read_id'].isin(fused_read_ids) == True]

    new_n50 = pd.concat([un_fused_df['sequence_length_template'], df2['combined_length']])
    stats_n50['New N50:'] = n50(new_n50)
    stats['Reads joined:'] = len(fused_read_ids)
    stats_n50['To be fused N50:'] = n50(split_df['sequence_length_template'])
    stats_n50['Fused read N50:'] = n50(df2['combined_length'])
    stats_n50['Un-fused reads N50:'] = n50(un_fused_df['sequence_length_template'])
    stats['Un-fused reads:'] = len(un_fused_df)
    stats['Fused reads:'] = len(df2)
    stats['New read count:'] = stats['Un-fused reads:'] + stats['Fused reads:']

    max_len = max([len(k) for k, v in stats.items()])
    for k, v in stats.items():
        print('{k:{m}}\t{v}'.format(k=k, m=max_len, v=v))
    # print("-------------------------------------------------")
    print("")
    stats_df = pd.DataFrame(stats_n50).T
    stats_df = stats_df[['MIN', 'MAX', 'MEAN', 'N50']]
    print(stats_df)
    print("")
    # print("-------------------------------------------------")
    # export fused reads summary file
    header = ['coords', 'filename', 'channel', 'start_time',
              'duration', 'combined_length', 'target_name', 'strand',
              'start_match', 'end_match', 'cat_read_id', 'count']
    df2.to_csv(args.out_fused, sep="\t", header=True, columns=header, index=False)
    print("Fused read summary file saved as {f}".format(f=args.out_fused))
    print("")
    print("Top ten fused reads by combined length:")
    top_n(df2, 'combined_length', 10)
    print("Top ten un-fused reads by length:")
    top_n(un_fused_df, 'sequence_length_template', 10)
    print("Top ten original reads by length:")
    top_n(ss, 'sequence_length_template', 10)
    """
    Output top 10 reads by length for both original and fused datasets
    fix start and end match for +/- strands
    make sure all vals are int/float accordingly
    'next_end', => 'end_time'
    'difference', => 'duration'
    'cat_read_id' => 
    """

    if args.debug:
        df2.to_csv('debug.csv', sep=",", index=False)


def top_n(df, field, n):
    df = df.sort_values(by=field, ascending=False)
    rows = df.filter([field], axis=1).head(n=n).reset_index()
    for idx, row in rows.iterrows():
        print('{i}:\t{len}'.format(i=idx+1, len=row[field]))
    return

def n50(lengths):
    """
    Calculate N50 from a list of lengths
    :param lengths: list of lengths as ints
    :return: N50 as an int
    """
    length_dict = {'length': lengths}
    df = pd.DataFrame(length_dict)
    df = df[df.length != 0]
    df = df.sort_values(by='length')
    df['cum_sum'] = df.length.cumsum()
    length_sum = df.length.sum()
    reads_n50 = int(df.length.where(df.cum_sum > length_sum / 2).dropna().iloc[0])
    reads_mean = int(df['length'].mean())
    reads_min = df['length'].min()
    reads_max = df['length'].max()
    res_dict = {'MIN': int(reads_min), 'MAX': int(reads_max), 'MEAN': int(reads_mean), 'N50': int(reads_n50)}
    return res_dict


def get_args():
    parser = ArgumentParser(
        description="""Parse sequencing_summary.txt files 
                       and .paf files to find split reads 
                       in an Oxford Nanopore Dataset""",
        add_help=False)
    general = parser.add_argument_group(
        title='General options')
    general.add_argument("-h", "--help",
                         action="help",
                         help="Show this help and exit"
                         )
    general.add_argument("-d", "--distance",
                         help='''Specify the maximum distance between consecutive mappings.
                              This is the difference between \'Target Start\' and \'Target End\' in 
                              the paf file. Defaults to 10000''',
                         type=int,
                         default=10000,
                         metavar=''
                         )
    general.add_argument("-D", "--debug",
                         help='''Write debug file''',
                         action="store_true",
                         default=False,
                         )
    in_args = parser.add_argument_group(
        title='Input sources'
    )
    in_args.add_argument("-s", "--summary",
                         help="A sequencing summary file generated by albacore",
                         type=str,
                         default="",
                         required=True,
                         metavar=''
                         )
    in_args.add_argument("-p", "--paf",
                         help="A paf file generated by minimap2",
                         type=str,
                         default='',
                         required=True,
                         metavar=''
                         )
    out_args = parser.add_argument_group(
        title='Output files'
    )
    out_args.add_argument('-F', '--out-fused',
                          help='''Specify name of the fused_read file. This file only contains chains of reads. 
                               Defaults to \'fused_reads.txt\'''',
                          type=str,
                          default='fused_reads.txt',
                          metavar=''
                          )
    return parser.parse_args()


def debug(args):
    dirs = dir(args)
    for attr in dirs:
        if attr[0] != '_':
            print('{a:<15} {b}'.format(a=attr, b=getattr(args, attr)))


if __name__ == '__main__':
    main()