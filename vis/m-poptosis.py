import pandas as pd
import numpy as np

input_arg = ""
df = pd.DataFrame.from_csv(input_arg, sep='\t', header=0)

#  1 - filename
#  2 - read_id
#  3 - run_id
#  4 - channel
#  5 - start_time
#  6 - duration
#  7 - num_events
#  8 - passes_filtering
#  9 - template_start
# 10 - num_events_template
# 11 - template_duration
# 12 - num_called_template
# 13 - sequence_length_template
# 14 - mean_qscore_template
# 15 - strand_score_template
# 16 - calibration_strand_genome_template
# 17 - calibration_strand_identity_template
# 18 - calibration_strand_accuracy_template
# 19 - aligned_speed_bps_template

df = df.sort_values(['run_id', 'channel', 'start_time'])

print(len(df))

# list of run ids
run_ids = df.run_id.drop_duplicates(keep="first")
run_id_dict = {}
for index, run_id in run_ids.items():
    # per run loop i.e. each restart
    run_id_dict[run_id] = df.where(df.run_id == run_id)
    run_id_dict[run_id] = run_id_dict[run_id].dropna()


plot_df = pd.DataFrame()

for key, value in run_id_dict.items():
    # create list of channels in each run_id
    channels = value.channel.drop_duplicates(keep="first")
    channels = channels.dropna()
    for index, ch in channels.items():
        # ndarray of values for ch in run_id
        not_last = value.sequence_length_template.where(value.channel == ch).dropna()[:-1].values
        last = value.sequence_length_template.where(value.channel == ch).dropna()[-1:].values
        plot_df = plot_df.append(pd.DataFrame({'length': not_last, 'not_last': np.ones(len(not_last), dtype=np.int)}))
        plot_df = plot_df.append(pd.DataFrame({'length': last, 'not_last': np.zeros(len(last), dtype=np.int)}))
        # get the list
        # channel = channel.sort_values(['start_time'])


print(len(plot_df))
