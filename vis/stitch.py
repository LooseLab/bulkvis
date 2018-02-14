import h5py
import pandas as pd

input_arg = "/Users/Alex/projects/bokehApp/bulkvis/data/NA12878_Run6_Sample2_29123.fast5"
output_arg = "/Users/Alex/projects/bokehApp/bulkvis/data/make/PLSP57501_20180118_FAH41614_MN18458_sequencing_run_NA12878_Run6_Sample2_29123_read_234_ch_391_strand.fast5"

bulkfile = h5py.File(input_arg, "r")
readfile = h5py.File(output_arg, "w")

version_num = 0.6

read_id_str = "b45a4b09-6f22-40f6-afd9-aa7fca8e89f3"
ch_num = 391
ch_str = "Channel_{ch}".format(ch=ch_num)

ugk = readfile.create_group("UniqueGlobalKey")

bulkfile.copy('UniqueGlobalKey/context_tags', ugk)
bulkfile.copy('UniqueGlobalKey/tracking_id', ugk)
bulkfile.copy("IntermediateData/{ch}/Meta".format(ch=ch_str), ugk)

readfile["UniqueGlobalKey"]["channel_id"] = readfile["UniqueGlobalKey"]["Meta"]
del readfile["UniqueGlobalKey"]["Meta"]

readfile["UniqueGlobalKey"]["channel_id"].attrs.create('channel_number', ch_num, None, dtype='<S4')
remove_attrs = ["description", "elimit", "scaling_used", "smallest_event", "threshold", "window"]
for attr in remove_attrs:
    del readfile["UniqueGlobalKey"]["channel_id"].attrs[attr]

int_data_path = bulkfile["IntermediateData"][ch_str]["Reads"]
int_dict = {
    'read_id': int_data_path["read_id"],
    'read_number': int_data_path["read_number"],
    'read_start': int_data_path["read_start"],
    'median_before': int_data_path["median_before"],
    'current_well_id': int_data_path["current_well_id"]
}
df = pd.DataFrame(data=int_dict)
df.read_id = df.read_id.str.decode('utf8')
df = df.where(df.read_id == read_id_str).dropna()
data_start = int(df.iloc[0].read_start)
data_end = int(df.iloc[-1].read_start) + 1737
read_number = int(df.iloc[0].read_number)
attrs = {
    'duration': {'val': data_end - data_start, 'd': 'uint32'},
    'median_before': {'val': df.iloc[0].median_before, 'd': 'Float64'},
    'read_id': {'val': read_id_str, 'd': '<S38'},
    'read_number': {'val': read_number, 'd': 'uint16'},
    'start_mux': {'val': int(df.iloc[0].current_well_id), 'd': 'uint8'},
    'start_time': {'val': data_start, 'd': 'uint64'}
}

dataset = bulkfile["Raw"][ch_str]["Signal"][()]
dataset = dataset[data_start:data_end]

readfile.create_group('Raw/Reads/Read_{n}'.format(n=read_number))
readfile.attrs.create('file_version', version_num, None, dtype='Float64')
# add read_### attrs
for k, v in attrs.items():
    readfile["Raw"]["Reads"]["Read_{n}".format(n=read_number)].attrs.create(k, v['val'], None, dtype=v['d'])

ms = [18446744073709551615]
ds = readfile.create_dataset('Raw/Reads/Read_{n}/Signal'.format(n=read_number), data=(dataset), maxshape=(ms), chunks=True, dtype='int16', compression="gzip", compression_opts=1)

bulkfile.close()
readfile.close()


"""
ToDo:
 - build 'Read_###' attributes
 - Fix dataset length issues
 - Configure plot options -> argparse
 - Errors...
"""