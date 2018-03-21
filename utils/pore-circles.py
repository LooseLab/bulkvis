from collections import defaultdict
from datetime import datetime

import channel_lookup
import h5py
import matplotlib.pyplot as plt
from bokeh.models import HoverTool, OpenURL, TapTool

from bokeh.io import curdoc
from bokeh.layouts import row, column, widgetbox
from bokeh.models.widgets import Select

from bokeh.plotting import figure, output_file, show, ColumnDataSource
from matplotlib.colors import rgb2hex

"""
Takes an input as a fast5 file and loops through the state data 
group to show/generate all the pores that have data recorded
"""


def group_counter(h5_file):
    group_count = 0
    if isinstance(h5_file, h5py.Group):
        for h5_group in h5_file:
            if isinstance(h5_file[h5_group], h5py.Group):
                group_count += 1
    return group_count


def enumerate_as_dict(path, field):
    """
    Enumerate a field datatype as a dictionary
    :param path: string in form 'fast5["path"]["to"]["dataset"]'
    :param field: string containing the field to get dtype for
    :return: dictionary of enumerated value and description
    """
    dataset_dtype = h5py.check_dtype(enum=path.dtype[field])
    # inv_dtype may lose some dataset dtypes there are duplicates of 'v'
    inv_dtype = {v: k for k, v in dataset_dtype.items()}
    return inv_dtype


def normalise_strand(n, min, max):
    norm_res = (n - min) / (max - min)
    return norm_res


input_arg = 'bulkfiles/PLSP61583_20171106_FAH20412_MN18458_sequencing_run_Notts_cDNA_Run1_96243.fast5'

file = h5py.File(input_arg, 'r')
time = datetime.utcnow()
time = time.strftime("%Y%m%d_%H%M%S")
fn = file.filename.split("/")[-1]
data = dict(
    x=[],
    y=[],
    channel=[],
    shape=[],
    state_below=[],
    state_inrange=[],
    state_above=[],
    state_good_single=[],
    state_strand=[],
    state_multiple=[],
    state_unavailable=[],
    state_unblocking=[],
    state_saturated=[],
    state_adapter=[],
    state_unclassified=[],
    state_unclassified_following_reset=[],
    state_pending_manual_reset=[],
    state_pending_mux_change=[],
    colour=[],
)

x_line = 1
y_line = 1
# channel_range = group_counter(file["Raw"]) + 1

# sample_frequency = int(file["UniqueGlobalKey"]["context_tags"].attrs["sample_frequency"].decode('utf8'))

statedata = file["StateData"]
pore = defaultdict()
channel_options = []
qi = 0
for group in statedata:
    channel_number = int(group.split("_")[1])
    # channel_options.append((group, channel_number))
    dataset = statedata[group]["States"]
    pore[group] = defaultdict(int)
    pore[group]["enum"] = enumerate_as_dict(dataset, "summary_state")
    if qi == 0:
        for k, v in pore[group]["enum"].items():
            channel_options.append(v)
        qi += 1
    for i in dataset['summary_state']:
        pore[group][i] += 1

    data['x'].append(channel_lookup.lookup[channel_number][0])
    data['y'].append(channel_lookup.lookup[channel_number][1])
    data['channel'].append(channel_number)
    data['shape'].append(dataset.len())
    data['state_below'].append(pore[group][0])
    data['state_inrange'].append(pore[group][1])
    data['state_above'].append(pore[group][2])
    data['state_good_single'].append(pore[group][3])
    data['state_strand'].append(pore[group][4])
    data['state_multiple'].append(pore[group][5])
    data['state_unavailable'].append(pore[group][6])
    data['state_unblocking'].append(pore[group][7])
    data['state_saturated'].append(pore[group][10])
    data['state_adapter'].append(pore[group][11])
    data['state_unclassified'].append(pore[group][200])
    data['state_unclassified_following_reset'].append(pore[group][201])
    data['state_pending_manual_reset'].append(pore[group][202])
    data['state_pending_mux_change'].append(pore[group][203])

# channel_options = sorted(channel_options, key=lambda x: x[1])
# channel_options = [i[0] for i in channel_options]
data_select = Select(title="Output:", value=channel_options[4], options=channel_options)


def colour_selector(colour):
    col_option = data[colour]
    strand_min = min(col_option)
    strand_max = max(col_option)
    colour = plt.get_cmap('Blues')

    for element in col_option:
        n_norm_res = normalise_strand(element, strand_min, strand_max)
        col_rgb = colour(n_norm_res)[:3]
        col_hex = rgb2hex(col_rgb)
        data['colour'].append(col_hex)

source = ColumnDataSource(data=data)
hover = HoverTool(tooltips=[
    ("Channel", "@channel"),
    ("Events", "@shape"),
    ("below", "@state_below"),
    ("inrange", "@state_inrange"),
    ("above", "@state_above"),
    ("good single", "@state_good_single"),
    ("strand", "@state_strand"),
    ("multiple", "@state_multiple"),
    ("unavailable", "@state_unavailable"),
    ("unblocking", "@state_unblocking"),
    ("saturated", "@state_saturated"),
    ("adapter", "@state_adapter"),
    ("unclassified", "@state_unclassified"),
    ("unclassified following reset", "@state_unclassified_following_reset"),
    ("pending manual reset", "@state_pending_manual_reset"),
    ("pending mux change", "@state_pending_mux_change"),
])

def update():
    select_value = "state_{}".format(data_select.value)
    cols = colour_selector(select_value)
    del source.data[colour][:]
    source.data = dict(colour=cols)


data_select.on_change('value', lambda attr, old, new: update())

p = figure(
    title="MinION Channels: " + fn,
    tools=[hover, "crosshair", "save", "tap"],
    # sizing_mode="stretch_both"
    plot_height=600,
    plot_width=800,
)

url = "#@channel"
taptool = p.select(type=TapTool)
taptool.callback = OpenURL(url=url)

p.text(
    data['x'],
    data['y'],
    text=data['channel'],
    text_font_size="10pt",
    text_baseline="middle",
    text_align="center"
)

r = p.circle(
    'x',
    'y',
    size=30,
    fill_color='colour',
    source=source
)

glyph = r.glyph
glyph.fill_alpha = 0.8
glyph.line_width = 1

p.axis.visible = False
p.grid.grid_line_color = None
p.toolbar.logo = None
p.toolbar_location = None
p.outline_line_width = 10
p.outline_line_alpha = 0.3
p.outline_line_color = "navy"

p.background_fill_color = "beige"
p.background_fill_alpha = 0.5

# output_file(time + "channel_diagram.html")

# show(p)
inputs = widgetbox(data_select)

try:
    update()
except Exception as e:
    print(e)

curdoc().add_root(row(inputs, p, sizing_mode="stretch_both"))
curdoc().title = "Channel Diagram"


file.close()
