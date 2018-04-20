from argparse import ArgumentParser
import numpy as np
import pandas as pd
import os.path
import sys
from whale_watch import fuse_reads


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
    df2, ss, chained_read_ids = fuse_reads(args.summary, args.paf, args.distance, 10, True, False)
    read_dict = {}

    #open file for writing
    file = open(args.out_fused, "w")
    myreadtracker = ReadTracker()
    ssclean = ss[ss['sequence_length_template'] != 0]
    print("{n} reads to process.".format(n=len(ssclean)))
    myreadtracker.settotal(len(ssclean))
    myreadtracker.settotal(len(ssclean))
    cnt = 1
    for dirpath, dirnames, filenames in os.walk(args.readfiles):
        for filename in [f for f in filenames if f.endswith(".fastq")]:
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
                            finished = 0
                            if len(df2['cat_read_id'].loc[df2['read_id'] == read_id].values[0].split('|')) > 1:
                                myreadtracker.fusedseen()
                                for value in df2['cat_read_id'].loc[df2['read_id'] == read_id].values[0].split('|'):
                                    if value in read_dict.keys():
                                        finished = 1
                                    else:
                                        finished = 0
                                        break
                                if finished == 1:
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
                                            readhead = "@{}".format(df2['cat_read_id'].loc[df2['read_id'] == value].values[0])
                                            readheadelements = read_dict[value]['header'].split()
                                            readheadelements[0] = readhead
                                            readhead = " ".join(readheadelements)
                                        read_dict.pop(value, None)
                                    myreadtracker.readwritten()
                                    file.write(readhead+"\n")
                                    file.write(readseq+"\n")
                                    file.write(readdiv+"\n")
                                    file.write(readqual+"\n")
                                else:
                                    pass
                        elif len(df2.loc[df2['cat_read_id'].str.contains(read_id)]) > 0:
                            finished = 0
                            if len(df2['cat_read_id'].loc[df2['cat_read_id'].str.contains(read_id)].values[0].split('|')) > 1:
                                myreadtracker.fusedseen()
                                for value in df2['cat_read_id'].loc[df2['cat_read_id'].str.contains(read_id)].values[0].split('|'):
                                    if value in read_dict.keys():
                                        finished = 1
                                    else:
                                        finished = 0
                                        break
                                if finished == 1:
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
                                            readhead = "@{}".format(df2['cat_read_id'].loc[df2['read_id'] == value].values[0])
                                            readheadelements = read_dict[value]['header'].split()
                                            readheadelements[0] = readhead
                                            readhead = " ".join(readheadelements)
                                        read_dict.pop(value, None)
                                    myreadtracker.readwritten()
                                    file.write(readhead+"\n")
                                    file.write(readseq+"\n")
                                    file.write(readdiv+"\n")
                                    file.write(readqual+"\n")
                                else:
                                    # Not found everything yet
                                    pass
                        else:
                            if read_id in chained_read_ids:
                                print("we should have found this")
                                sys.exit()
                            myreadtracker.unfusedseen()
                            if args.W:
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


def get_args():
    parser = ArgumentParser(
        description="""Parse sequencing_summary.txt files and .paf files to find chained reads 
                       in an Oxford Nanopore Dataset and output fused fastq files""",
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
                         help="Full path to the folder containing fastq files you wish to join",
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
    out_args.add_argument('-W',
                          help='''Outputs just the fused reads''',
                          action="store_false",
                          default=True,
                          metavar=''
                          )
    return parser.parse_args()


if __name__ == '__main__':
    main()
