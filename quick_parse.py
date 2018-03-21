import numpy as np
import pandas as pd


def n50(lengths):
    """
    Calculate N50 from a list of lengths
    :param lengths: list of lengths as ints
    :return: N50 as an int
    """
    length_dict = {'length': lengths}
    df = pd.DataFrame(length_dict)
    df = df.sort_values(by='length')
    df['cum_sum'] = df.length.cumsum()
    length_sum = df.length.sum()
    n50 = int(df.length.where(df.cum_sum > length_sum / 2).dropna().iloc[0])
    mean = int(df['length'].mean())
    min = df['length'].min()
    max = df['length'].max()
    return min, max, mean, n50


sequencing_summary = "sequencing_summary.txt"
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

old_n50 = n50(ss['sequence_length_template'])
print("Original N50:\t{v}".format(v=old_n50))

paf_file = "allreads.paf"
pf = pd.read_csv(paf_file, sep='\t', header=None, usecols=[0,4,5,7,8])
pf.columns = ['Qname','Strand','Tname','Tstart','Tend']
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
df2 = df2[df2['match_distance'] < 10000]
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

df2['start_match']=df2.groupby(level='cs')['Tstart_A'].first()
df2['end_match']=df2.groupby(level='cs')['Tend_B'].last()

df2['start_time'] = df2.groupby(level='cs')['start_time'].first()
df2['next_end'] = df2.groupby(level='cs')['next_end'].last()
df2['difference'] = df2['next_end']-df2['start_time']
df2['stime_floor'] = np.floor(df2['start_time']).astype('int64').astype('str')
df2['etime_ceil'] = np.ceil(df2['next_end']).astype('int64').astype('str')
df2['coords'] = df2['channel'].astype('int64').astype('str') + ":" + df2['stime_floor'] + "-" + df2['etime_ceil']

df2['filename'] = df2['filename'].str.split('_read', 1).str[0]
df2['combined_length'] = df2['combined_length'] + df2['last_length']

df2['cat_read_id'] = df2['cat_read_id'] + "|" + df2['last_read_id']
header = ['coords', 'filename', 'channel', 'start_time', 'next_end','difference', 'combined_length','Tname_A','Strand_A', 'start_match', 'end_match','cat_read_id']

chained_read_ids = pd.concat([df2['read_id'], df2['next_read_id']])
df2 = df2.drop_duplicates(subset=header, keep='first')

chained_read_ids = chained_read_ids.unique()


ss2 = ss[ss['read_id'].isin(chained_read_ids) == False]

ss3 = ss[ss['read_id'].isin(chained_read_ids) == True]

reads_to_be_fused_N50 = n50(ss3['sequence_length_template'])
unjoined_reads_N50 = n50(ss2['sequence_length_template'])
fused_reads_N50 = n50(df2['combined_length'])

new_n50 = pd.concat([ss2['sequence_length_template'], df2['combined_length']])
new_n50_num = n50(new_n50)


print("New N50:\t{v}".format(v=new_n50_num))



print ("No. reads joined:\t{v}".format(v=len(chained_read_ids)))
print ("To be fused N50:\t{v}".format(v=reads_to_be_fused_N50))
print ("Fused read N50:\t{v}".format(v=fused_reads_N50))
print ("unfused reads N50:\t{v}".format(v=unjoined_reads_N50))
print ("Fused reads:\t{v}".format(v=len(df2)))

# Begin export here
# df2.to_csv('complete_df.csv', sep=",", index=False)
df2.to_csv('fused_reads_summary.txt', sep="\t", header=True, columns=header, index=False)
