from argparse import ArgumentParser
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import sys


class ReadTracker():
    def __init__(self):
        #print("ReadTracker initiated")
        self.readdict=dict()
        self.readdict["total"]=0
        self.readdict["fused"]=0
        self.readdict["unfused"]=0
        self.readdict["final_total"]=0

    def settotal(self,total):
        self.readdict["total"]=total

    def fusedseen(self):
        self.readdict["fused"]+=1
        self.readseen()

    def unfusedseen(self):
        self.readdict["unfused"]+=1
        self.readseen()
        self.readwritten()

    def readseen(self):
        self.readdict["total"]-=1

    def readwritten(self):
        self.readdict["final_total"]+=1

    def result(self):
        return "{} to process, {} fused, {} unfused, {} written.".format(self.readdict["total"],self.readdict["fused"],self.readdict["unfused"],self.readdict["final_total"])


def main():
    args = get_args()

    sequencing_summary = args.summary
    ss_fields = ['channel', 'start_time', 'duration', 'run_id', 'read_id', 'sequence_length_template', 'filename']
    ss = pd.read_csv(sequencing_summary, sep='\t', usecols=ss_fields)
    ss = ss.sort_values(by=['channel', 'run_id', 'start_time'])
    ss['diff'] = ss['start_time'].shift(-1) - (ss['start_time'] + ss['duration'])
    # quickly put next read_id on as another column
    ss['next_read_id'] = ss['read_id'].shift(-1)
    ss['next_start_time'] = ss['start_time'].shift(-1)
    ss['next_end'] = ss['next_start_time'] + ss['duration'].shift(-1)
    ss['next_sequence_length_template'] = ss['sequence_length_template'].shift(-1)
    ss['combined_length'] = ss['sequence_length_template'] + ss['next_sequence_length_template']
    ss['next_sequence_length_template'] = ss['next_sequence_length_template'].fillna(0).astype('int64')
    ss['combined_length'] = ss['combined_length'].fillna(0).astype('int64')

    stats = OrderedDict()
    stats['Original N50:'] = n50(ss['sequence_length_template'])

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

    cond_1 = df2['next_read_id'] == df2['read_id'].shift(-1)
    cond_2 = df2['read_id'] == df2['next_read_id'].shift(-1)

    df2['COND'] = np.where(cond_1 | cond_2, True, False)
    df2['W'] = np.where(df2['COND'].shift(1) == False, 1, 0)
    df2['cs'] = df2['W'].cumsum()
    df2 = df2.set_index('cs')

    df2['cat_read_id'] = df2.groupby(level='cs')['read_id'].apply('|'.join)
    df2['last_read_id'] = df2.groupby(level='cs')['next_read_id'].last()

    df2['combined_length'] = df2.groupby(level='cs')['sequence_length_template'].sum()
    df2['last_length'] = df2.groupby(level='cs')['next_sequence_length_template'].last()

    df2['start_match'] = df2.groupby(level='cs')['Tstart_A'].first()
    df2['end_match'] = df2.groupby(level='cs')['Tend_B'].last()

    df2['start_time'] = df2.groupby(level='cs')['start_time'].first()
    df2['next_end'] = df2.groupby(level='cs')['next_end'].last()
    df2['difference'] = df2['next_end'] - df2['start_time']
    df2['stime_floor'] = np.floor(df2['start_time']).astype('int64').astype('str')
    df2['etime_ceil'] = np.ceil(df2['next_end']).astype('int64').astype('str')
    df2['coords'] = df2['channel'].astype('int64').astype('str') + ":" + df2['stime_floor'] + "-" + df2['etime_ceil']

    df2['filename'] = df2['filename'].str.split('_read', 1).str[0]
    df2['combined_length'] = df2['combined_length'] + df2['last_length']

    df2['cat_read_id'] = df2['cat_read_id'] + "|" + df2['last_read_id']

    header = ['coords', 'channel', 'start_time', 'next_end',
              'difference', 'combined_length', 'start_match',
              'end_match', 'cat_read_id']

    chained_read_ids = pd.concat([df2['read_id'], df2['next_read_id']])

    df2 = df2.drop_duplicates(subset=header, keep='first')

    chained_read_ids = chained_read_ids.unique()
    ss2 = ss[ss['read_id'].isin(chained_read_ids) == False]
    ss3 = ss[ss['read_id'].isin(chained_read_ids) == True]

    new_n50 = pd.concat([ss2['sequence_length_template'], df2['combined_length']])
    stats['New N50:'] = n50(new_n50)
    stats['Reads joined:'] = len(chained_read_ids)
    stats['To be fused N50:'] = n50(ss3['sequence_length_template'])
    stats['Fused read N50:'] = n50(df2['combined_length'])
    stats['Un-fused reads N50:'] = n50(ss2['sequence_length_template'])
    stats['Fused reads:'] = len(df2)


    max_len = max([len(k) for k, v in stats.items()])
    for k, v in stats.items():
        print('{k:{m}}\t{v}'.format(k=k, m=max_len, v=v))

    #df2.to_csv(args.out_summary, sep="\t", header=True, columns=header, index=False)
    read_dict = {}

    #open file for writing
    file = open(args.out_fused, "w")

    myreadtracker = ReadTracker()

    ssclean = ss[ss['sequence_length_template'] != 0]
    print("{n} reads to process.".format(n=len(ssclean)))
    myreadtracker.settotal(len(ssclean))
    myreadtracker.settotal(len(ssclean))

    #check reads
    cnt = 1
    for dirpath, dirnames, filenames in os.walk(args.readfiles):
        for filename in [f for f in filenames if f.endswith(".fastq")]:
            #print(os.path.join(dirpath, filename))
            filepath = os.path.join(dirpath, filename)
            with open(filepath) as fp:
                line = fp.readline()

                while line:
                    if len(line.split()) > 1: #means we have a fastq ID line
                        if cnt%400 == 0:
                            print(myreadtracker.result())
                        read_id = line.split()[0][1:]

                        read_dict[read_id] = dict()
                        read_dict[read_id]['header'] = line.strip()
                        read_dict[read_id]['fasta'] = fp.readline().strip()
                        read_dict[read_id]['divide'] = fp.readline().strip()
                        read_dict[read_id]['quality'] = fp.readline().strip()
                        if len(df2.loc[df2['read_id'] == read_id]) > 0:
                            #print(df2.loc[df2['read_id'] == read_id])
                            #print("Line {}: {}".format(cnt, line.strip()))
                            finished = 0
                            if len(df2['cat_read_id'].loc[df2['read_id'] == read_id].values[0].split('|')) > 1:
                                #print(len(df2['cat_read_id'].loc[df2['read_id'] == read_id].values[0].split('|')))
                                myreadtracker.fusedseen()
                                for value in df2['cat_read_id'].loc[df2['read_id'] == read_id].values[0].split('|'):
                                    #print(value)
                                    if value in read_dict.keys():
                                        finished = 1
                                        #print(finished)
                                    else:
                                        finished = 0
                                        break
                                if finished == 1:
                                    #print("We have all we need!")
                                    readhead = ""
                                    readseq = ""
                                    readdiv = ""
                                    readqual = ""
                                    for value in df2['cat_read_id'].loc[df2['read_id'] == read_id].values[0].split('|'):
                                        readseq = readseq + read_dict[value]['fasta']
                                        readqual = readqual + read_dict[value]['quality']
                                        if len(readdiv) == 0:
                                            readdiv = read_dict[value]['divide']
                                        if len(readhead) == 0:
                                            #@5e8cd22a-d788-42ad-a527-3ee7ef244f23 runid=df66bb2ee79201ee2b42723347c86856f0204b6d read=1676 ch=243 start_time=2018-01-27T02:08:35Z
                                            readhead = "@{}".format(df2['cat_read_id'].loc[df2['read_id'] == value].values[0])
                                            readheadelements = read_dict[value]['header'].split()
                                            readheadelements[0] = readhead
                                            readhead = " ".join(readheadelements)
                                        #print(readhead)
                                        read_dict.pop(value, None)
                                    myreadtracker.readwritten()
                                    file.write(readhead+"\n")
                                    file.write(readseq+"\n")
                                    file.write(readdiv+"\n")
                                    file.write(readqual+"\n")
                                else:
                                    #print("Not found everything yet")
                                    pass
                        elif len(df2.loc[df2['cat_read_id'].str.contains(read_id)]) > 0:
                            #print("this is where we should have found it - we had {} reads".format(len(df2.loc[df2['cat_read_id'].str.contains(read_id)])))
                            finished = 0
                            if len(df2['cat_read_id'].loc[df2['cat_read_id'].str.contains(read_id)].values[0].split('|')) > 1:
                                # print(len(df2['cat_read_id'].loc[df2['read_id'] == read_id].values[0].split('|')))
                                myreadtracker.fusedseen()
                                for value in df2['cat_read_id'].loc[df2['cat_read_id'].str.contains(read_id)].values[0].split('|'):
                                    #print(value)
                                    if value in read_dict.keys():
                                        finished = 1
                                        #print(finished)
                                    else:
                                        finished = 0
                                        break
                                if finished == 1:
                                    #print("We have all we need!")
                                    readhead = ""
                                    readseq = ""
                                    readdiv = ""
                                    readqual = ""
                                    for value in df2['cat_read_id'].loc[df2['cat_read_id'].str.contains(read_id)].values[0].split('|'):
                                        readseq = readseq + read_dict[value]['fasta']
                                        readqual = readqual + read_dict[value]['quality']
                                        if len(readdiv) == 0:
                                            readdiv = read_dict[value]['divide']
                                        if len(readhead) == 0:
                                            # @5e8cd22a-d788-42ad-a527-3ee7ef244f23 runid=df66bb2ee79201ee2b42723347c86856f0204b6d read=1676 ch=243 start_time=2018-01-27T02:08:35Z
                                            readhead = "@{}".format(df2['cat_read_id'].loc[df2['read_id'] == value].values[0])
                                            readheadelements = read_dict[value]['header'].split()
                                            readheadelements[0] = readhead
                                            readhead = " ".join(readheadelements)
                                        #print(readhead)
                                        read_dict.pop(value, None)
                                    myreadtracker.readwritten()
                                    file.write(readhead+"\n")
                                    file.write(readseq+"\n")
                                    file.write(readdiv+"\n")
                                    file.write(readqual+"\n")
                                    # sys.exit()
                                else:
                                    #print("Not found everything yet")
                                    pass
                            #sys.exit()
                            #df[df['ids'].str.contains("ball")]

                        else:
                            if read_id in chained_read_ids:
                                print("we should have found this")
                                sys.exit()
                            #print("This read is not a merged read.")
                            myreadtracker.unfusedseen()
                            file.write(read_dict[read_id]['header']+"\n")
                            file.write(read_dict[read_id]['fasta']+"\n")
                            file.write(read_dict[read_id]['divide']+"\n")
                            file.write(read_dict[read_id]['quality']+"\n")
                            read_dict.pop(read_id,None)
                            pass
                    line = fp.readline()
                    cnt += 1
    print("!!!!!!!", len(read_dict))
    print(myreadtracker.result())
    file.close()



def n50(lengths):
    """
    Calculate N50 from a list of lengths
    :param lengths: list of lengths as ints
    :return: N50 as an int
    """
    length_dict = {'length': lengths}
    df = pd.DataFrame(length_dict)
    df = df[df.length!=0]
    df = df.sort_values(by='length')
    df['cum_sum'] = df.length.cumsum()
    length_sum = df.length.sum()
    reads_n50 = int(df.length.where(df.cum_sum > length_sum / 2).dropna().iloc[0])
    reads_mean = int(df['length'].mean())
    reads_min = df['length'].min()
    reads_max = df['length'].max()
    return int(reads_min), int(reads_max), int(reads_mean), int(reads_n50)


def get_args():
    parser = ArgumentParser(
        description="""Parse sequencing_summary.txt files 
                       and .paf files to find chained reads 
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
    in_args.add_argument("-f", "--readfiles",
                         help="Full path to the folder containing fastq files you wish to join.",
                         type=str,
                         default='',
                         required=True,
                         metavar='')
    out_args = parser.add_argument_group(
        title='Output files'
    )
    out_args.add_argument('-o', '--out-fused',
                          help='''Specify name of the fused_read fastq file. This file will contain fused reads and the remaining singleton reads. 
                               Defaults to \'fused_reads.fastq\'''',
                          type=str,
                          default='fused_reads.fastq',
                          metavar=''
                          )
    return parser.parse_args()


if __name__ == '__main__':
    main()
