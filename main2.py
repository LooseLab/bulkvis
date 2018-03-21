"""
/Users/Alex/projects/bokehApp/bulkvis/data/NA12878_Run6_Sample2_29123.fast5

@b45a4b09-6f22-40f6-afd9-aa7fca8e89f3 runid=f9291b45b0c66faa77755e51738d193fcfafffc7 read=234 ch=391 start_time=2018-01-18T21:59:40Z
@b45a4b09-6f22-40f6-afd9-aa7fca8e89f3, 391
"""
from collections import OrderedDict
import configparser
import math
import os

import h5py
import numpy as np
import pandas as pd
from bokeh.plotting import curdoc, figure
from bokeh.layouts import row, widgetbox
from bokeh.models import TextInput, Button, Toggle, Div, Range1d, Label, Span, CheckboxGroup, Dropdown, PreText

config = configparser.ConfigParser()
config.read(os.path.dirname(os.path.realpath(__file__)) + '/config.ini')
cfg = config['plot_opts']


def update_file():
    """"""
    if app_data['bulkfile']:
        app_data['bulkfile'].flush()
        app_data['bulkfile'].close()

    file_src = app_data['wdg_dict']['data_input'].value
    app_data.clear()
    app_data['app_vars'] = {}
    app_data['wdg_dict'] = OrderedDict()
    app_data['label_dt'] = OrderedDict()
    app_data['file_src'] = file_src
    app_data['INIT'] = True

    (app_data['bulkfile'],
     app_data['app_vars']['sf'],
     app_data['app_vars']['channel_list']) = open_bulkfile(app_data['file_src'])

    raw_path = app_data['bulkfile']["Raw"]
    for i, member in enumerate(raw_path):
        if i == 0:
            signal_ds = raw_path[member]["Signal"][()]
            app_data['app_vars']['len_ds'] = math.ceil(len(signal_ds) / app_data['app_vars']['sf'])

    app_data['wdg_dict'] = init_wdg_dict()
    app_data['wdg_dict']['data_input'].value = file_src
    app_data['wdg_dict']['fastq_input'] = TextInput(title="FASTQ sequence ID:", value="", css_classes=[])
    app_data['wdg_dict']['channel_select'] = TextInput(title="Channel:", value="", css_classes=[])

    app_data['wdg_dict']['squiggle_start'] = TextInput(title='Start time (seconds):', value=str(0), css_classes=[])
    app_data['wdg_dict']['squiggle_duration'] = TextInput(title='Duration (seconds):', value=str(0), css_classes=[])

    app_data['wdg_dict']['fastq_input'].on_change("value", go_to_fastq)
    app_data['wdg_dict']['channel_select'].on_change("value", update)

    layout.children[0] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg['wdg_width']))


def open_bulkfile(path):
    """"""
    # !!! add in check to see if this is a ONT bulkfile
    file = h5py.File(path, "r")
    sf = int(file["UniqueGlobalKey"]["context_tags"].attrs["sample_frequency"].decode('utf8'))
    # make channel_list
    channel_list = np.arange(1, len(file["Raw"]) + 1, 1).tolist()
    return file, sf, channel_list


def build_widgets():
    """"""
    check_labels = []
    jump_list = []
    check_active = []
    app_data['label_mp'] = {}
    for k, v in enumerate(app_data['label_dt'].items()):
        app_data['label_mp'][v[0]] = k
        check_labels.append(v[1])
        check_active.append(k)
        jump_list.append((v[1], str(v[0])))


    wdg = app_data['wdg_dict']
    wdg['channel_select'] = TextInput(title='Channel:', value=str(app_data['app_vars']['channel_num']), css_classes=[])
    wdg['squiggle_start'] = TextInput(title='Start time (seconds):', value=str(0), css_classes=[])
    wdg['squiggle_duration'] = TextInput(title='Duration (seconds):', value=str(0), css_classes=[])

    wdg['jump_next'] = Dropdown(label="Jump to next", button_type="primary", menu=jump_list)
    wdg['jump_next_err'] = PreText(text="")
    wdg['jump_prev'] = Dropdown(label="Jump to previous", button_type="primary", menu=jump_list)

    wdg['toggle_x_axis'] = Toggle(
        label="Fixed x-axis",
        button_type="danger",
        css_classes=['toggle_button_g_r'],
        active=True
    )
    wdg['toggle_y_axis'] = Toggle(
        label="Fixed Y-axis",
        button_type="danger",
        css_classes=['toggle_button_g_r'],
        active=True
    )
    wdg['toggle_annotations'] = Toggle(
        label="Display annotations",
        button_type="danger",
        css_classes=['toggle_button_g_r'],
        active=True
    )

    wdg['label_options'] = Div(text='Select annotations', css_classes=['filter-dropdown', 'caret-down'])
    wdg['label_filter'] = CheckboxGroup(labels=check_labels, active=check_active, css_classes=['filter-drop'])

    wdg['plot_options'] = Div(text='Plot Adjustments', css_classes=['adjust-dropdown', 'caret-down'])
    wdg['po_width'] = TextInput(title='Plot Width (px)', value=cfg['plot_width'], css_classes=['adjust-drop'])
    wdg['po_height'] = TextInput(title='Plot Height (px)', value=cfg['plot_height'], css_classes=['adjust-drop'])
    wdg['po_x_width'] = TextInput(title="x width", value=cfg['x_width'], css_classes=['adjust-drop'])
    wdg['po_y_max'] = TextInput(title="y max", value=cfg['y_max'], css_classes=['adjust-drop'])
    wdg['po_y_min'] = TextInput(title="y min", value=cfg['y_min'], css_classes=['adjust-drop'])
    wdg['label_height'] = TextInput(
        title="Annotation height (y-axis)",
        value=cfg['label_height'],
        css_classes=['adjust-drop']
    )
    wdg['toggle_smoothing'] = Toggle(
        label="Toggle smoothing",
        button_type="danger",
        css_classes=['toggle_button_o_r', 'adjust-drop'],
        active=True
    )

    wdg['channel_select'].on_change('value', update)
    wdg['label_filter'].on_change('active', update)
    wdg['jump_next'].on_click(next_update)
    wdg['jump_prev'].on_click(prev_update)

    for name in toggle_inputs:
        wdg[name].on_click(toggle_button)
    for name in int_inputs:
        wdg[name].on_change('value', is_input_int)

    return wdg


def create_figure(x_data, y_data, label_df, label_dt, wdg, app_vars):
    if wdg["toggle_smoothing"].active:
        w_range = int(wdg['squiggle_duration'].value)
        thin_factor = thinning_factor(w_range)
    else:
        thin_factor = 4
    if thin_factor == 0:
        thin_factor = 1
    thin_factor = int(thin_factor)

    greater_delete_index = np.argwhere(y_data > 2200)
    x_data = np.delete(x_data, greater_delete_index)
    y_data = np.delete(y_data, greater_delete_index)

    lesser_delete_index = np.argwhere(y_data < -1000)
    x_data = np.delete(x_data, lesser_delete_index)
    y_data = np.delete(y_data, lesser_delete_index)

    p = figure(
        plot_height=int(wdg['po_height'].value),
        plot_width=int(wdg['po_width'].value),
        title="{ch} Raw Output at {sf} events per second: {tf}".format(
            ch=app_vars['channel_str'],
            sf=app_vars['sf'],
            tf=thin_factor
        ),
        toolbar_location="above",
        tools=['xpan', 'xbox_zoom', 'save', 'reset'],
    )
    p.yaxis.axis_label = "Current (pA)"
    p.yaxis.major_label_orientation = "horizontal"
    p.xaxis.axis_label = "Time (seconds)"
    p.line(x_data[::thin_factor], y_data[::thin_factor], line_width=1)
    if wdg['toggle_y_axis'].active:
        p.y_range = Range1d(int(wdg['po_y_min'].value), int(wdg['po_y_max'].value))

    if wdg['toggle_x_axis'].active:
        p.x_range = Range1d(
            int(int(wdg['squiggle_start'].value)),
            int(wdg['squiggle_start'].value) + int(wdg['po_x_width'].value)
        )

    p.xaxis.major_label_orientation = math.radians(45)

    if wdg['toggle_annotations'].active:
        slim_label_df = label_df[(label_df['read_start'] >= int(wdg['squiggle_start'].value)) &
                                 (label_df['read_start'] <= app_vars['end_time'])]

        for index, label in slim_label_df.iterrows():
            if app_data['label_mp'][label.modal_classification] in wdg['label_filter'].active:
                event_line = Span(
                    location=label.read_start,
                    dimension='height',
                    line_color='green',
                    line_dash='dashed',
                    line_width=1
                )
                p.add_layout(event_line)
                labels = Label(
                    x=label.read_start,
                    y=int(wdg['label_height'].value),
                    text=str(label_dt[label.modal_classification]),
                    level='glyph',
                    x_offset=0,
                    y_offset=0,
                    render_mode='canvas',
                    angle=-300
                )
                p.add_layout(labels)
    return p


def thinning_factor(window_range):
    divisor = math.e ** 2.5
    factor = window_range / divisor
    return math.ceil(factor)


def init_wdg_dict():
    wdg_dict = OrderedDict()
    wdg_dict['data_input'] = TextInput(
        title='Data Source:',
        value="/Users/Alex/projects/bokehApp/bulkvis/data/NA12878_Run6_Sample2_29123.fast5"
    )
    wdg_dict['data_button'] = Button(label="Update source file", button_type="primary")
    wdg_dict['data_button'].on_click(update_file)
    return wdg_dict


def update_data(start_time, sf, bulkfile, channel):
    duration = start_time[1] - start_time[0]
    # get times and squiggles
    end_time, start_squiggle, end_squiggle = get_time_frames(start_time[0], duration, sf)
    # get data in numpy arrays
    x_data = points_index(start_time[0], end_time, sf)
    y_data, len_ds = signal_slice(start_squiggle, end_squiggle, bulkfile, channel)
    len_ds = len_ds / sf
    # get annotations
    label_df, label_dtypes = get_annotations(bulkfile, channel, sf)
    return x_data, y_data, label_df, label_dtypes, end_time, start_squiggle, end_squiggle, len_ds


def get_time_frames(data_start, duration, sf):
    start_squiggle = math.floor(data_start * sf)
    end_time = data_start + duration
    end_squiggle = math.ceil(start_squiggle + sf * duration)
    return end_time, start_squiggle, end_squiggle


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


def points_index(start_time, end_time, sf):
    step = 1 / sf
    points_arr = np.arange(start_time, end_time, step)
    return points_arr


def signal_slice(data_start, data_end, file, channel):
    dataset = file["Raw"][channel]["Signal"][()]
    len_ds = len(dataset)
    return dataset[data_start:data_end], len_ds


def get_annotations(file, channel, sf):
    int_data_path = file["IntermediateData"][channel]["Reads"]
    int_data_dtypes = enumerate_as_dict(int_data_path, "modal_classification")
    int_data_labels = {
        'read_id': int_data_path["read_id"],
        'read_start': int_data_path["read_start"],
        'modal_classification': int_data_path["modal_classification"]
    }
    labels_df = pd.DataFrame(data=int_data_labels).drop_duplicates(subset="read_id", keep="first")
    labels_df.read_start = labels_df.read_start / sf
    labels_df.read_id = labels_df.read_id.str.decode('utf8')
    return labels_df, int_data_dtypes


def is_input_int(attr, old, new):
    try:
        int(new)
        for wdg in int_inputs:
            if (app_data['wdg_dict'][wdg].value == new) and ('input-error' in app_data['wdg_dict'][wdg].css_classes):
                input_error(app_data['wdg_dict'][wdg], 'remove')
    except ValueError:
        for wdg in int_inputs:
            if app_data['wdg_dict'][wdg].value == new:
                input_error(app_data['wdg_dict'][wdg], 'add')
                return

    new = new.lstrip('0')
    update(attr, old, new)


def update(attr, old, new):
    try:
        channel = int(app_data['wdg_dict']['channel_select'].value)
        if channel in app_data['app_vars']['channel_list']:
            input_error(app_data['wdg_dict']['channel_select'], 'remove')
        else:
            input_error(app_data['wdg_dict']['channel_select'], 'add')
            return
    except ValueError:
        input_error(app_data['wdg_dict']['channel_select'], 'add')
        return

    app_data['app_vars']['channel_num'] = channel
    app_data['app_vars']['channel_str'] = "Channel_{ch}".format(ch=channel)
    end_time = int(app_data['wdg_dict']['squiggle_start'].value) + int(app_data['wdg_dict']['squiggle_duration'].value)
    # app_data['wdg_dict'][''].value =
    (app_data['x_data'],
     app_data['y_data'],
     app_data['label_df'],
     app_data['label_dt'],
     app_data['app_vars']['end_time'],
     app_data['app_vars']['start_squiggle'],
     app_data['app_vars']['end_squiggle'],
     app_data['app_vars']['len_ds']) = update_data(
        (int(app_data['wdg_dict']['squiggle_start'].value), end_time),
         app_data['app_vars']['sf'],
         app_data['bulkfile'],
         app_data['app_vars']['channel_str']
    )
    if app_data['INIT']:
        build_widgets()
        layout.children[0] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg['wdg_width']))
        app_data['INIT'] = False

    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data['label_df'],
        app_data['label_dt'],
        app_data['wdg_dict'],
        app_data['app_vars']

    )


def go_to_fastq(attr, old, new):
    fq = app_data['wdg_dict']['fastq_input'].value
    if fq[0] == "@":
        fq_list = fq[1:]
        fq_list = fq_list.split(" ")
        read_id = fq_list[0]
        channel_number = fq_list[3].split("=")[1]
        input_error(app_data['wdg_dict']['fastq_input'], 'remove')
    else:
        input_error(app_data['wdg_dict']['fastq_input'], 'add')
        print("Not FASTQ string!")
        return

    app_data['app_vars']['channel_num'] = channel_number
    app_data['app_vars']['channel_str'] = "Channel_{ch}".format(ch=channel_number)
    int_data_path = app_data['bulkfile']["IntermediateData"][app_data['app_vars']['channel_str']]["Reads"]
    int_data_labels = {
        'read_id': int_data_path["read_id"],
        'read_start': int_data_path["read_start"],
    }
    df = pd.DataFrame(data=int_data_labels)
    df.read_start = df.read_start / app_data['app_vars']['sf']
    df.read_id = df.read_id.str.decode('utf8')
    df = df.where(df.read_id == read_id)
    df = df.dropna()
    # !!! check that multiple rows are still here
    fq_start_time = math.floor(df.iloc[0, :].read_start)
    fq_end_time = math.ceil(df.iloc[-1, :].read_start)

    # build widgets
    (app_data['x_data'],
     app_data['y_data'],
     app_data['label_df'],
     app_data['label_dt'],
     app_data['app_vars']['end_time'],
     app_data['app_vars']['start_squiggle'],
     app_data['app_vars']['end_squiggle'],
     app_data['app_vars']['len_ds']) = update_data(
        (fq_start_time, fq_end_time),
        app_data['app_vars']['sf'],
        app_data['bulkfile'],
        app_data['app_vars']['channel_str']
    )
    app_data['wdg_dict'] = build_widgets()
    app_data['wdg_dict']['squiggle_start'].value = str(fq_start_time)
    app_data['wdg_dict']['squiggle_duration'].value = str(fq_end_time - fq_start_time)
    layout.children[0] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg['wdg_width']))
    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data['label_df'],
        app_data['label_dt'],
        app_data['wdg_dict'],
        app_data['app_vars']
    )


def toggle_button(state):
    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data['label_df'],
        app_data['label_dt'],
        app_data['wdg_dict'],
        app_data['app_vars']
    )


def input_error(widget, mode):
    """"""
    if mode == 'add':
        widget.css_classes.append('input-error')
    elif mode == 'remove':
        if widget.css_classes:
            del widget.css_classes[-1]
    else:
        print("mode not recognised")


def next_update(value):
    value = int(value)
    jump_start = app_data['label_df'][(app_data['label_df']['read_start'] > int(app_data['wdg_dict']['squiggle_start'].value) + 1) &
                                      (app_data['label_df']['modal_classification'] == value)]
    app_data['wdg_dict']['squiggle_start'].value = str(math.floor(jump_start['read_start'].iloc[0]))
    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data['label_df'],
        app_data['label_dt'],
        app_data['wdg_dict'],
        app_data['app_vars']
    )
    app_data['wdg_dict']['jump_next'].value = ""


def prev_update(value):
    value = int(value)
    jump_start = app_data['label_df'][(app_data['label_df']['read_start'] < int(app_data['wdg_dict']['squiggle_start'].value)) &
                                      (app_data['label_df']['modal_classification'] == value)]
    app_data['wdg_dict']['squiggle_start'].value = str(math.floor(jump_start['read_start'].iloc[-1]))
    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data['label_df'],
        app_data['label_dt'],
        app_data['wdg_dict'],
        app_data['app_vars']
    )
    app_data['wdg_dict']['jump_prev'].value = ""


def parse_coordinates(position):
    """
    Parse coordinate string into channel, start and end
    :param position: string, "channel:start-end" eg: "391:200-800"
    :return: str, int, int
    """
    coords = position.split(":")
    times = coords[1].split("-")
    channel_str = "Channel_{num}".format(num=coords[0])
    (start, end) = int(times[0]), int(times[1])
    return channel_str, start, end


app_data = {
    'file_src': None,               # bulkfile path (string)
    'bulkfile': None,               # bulkfile object
    'x_data': None,                 # numpy ndarray time points
    'y_data': None,                 # numpy ndarray signal data
    'label_df': None,               # pandas df of signal labels
    'label_dt': None,               # dict of signal enumeration
    'label_mp': None,               # dict matching labels to widget filter
    'app_vars': {                   # dict of variables used in plots and widgets
        'len_ds': None,                 # length of signal dataset in (seconds or events)
        'start_time': None,             # squiggle start time in seconds
        'end_time': None,               # squiggle end time in seconds
        'start_squiggle': None,         # squiggle start position (events)
        'end_squiggle': None,           # squiggle end position (events)
        'channel_str': None,            # 'Channel_NNN' (string)
        'channel_num': None,            # Channel number (int)
        'duration': None,               # squiggle duration (seconds)
        'sf': None,                     # sample frequency (int)
        'channel_list': None,           # list of all channels as int
        'fastq_str': None,              # fastq sequence id (string)
        'fastq_read_id': None           # fastq read id (string)
    },
    'wdg_dict': None,               # dictionary of widgets
    'controls': None,               # widgets added to widgetbox
    'pore_plt': None,               # the squiggle plot
    'INIT': True                    # Initial plot with bulkfile (bool)
}

int_inputs = ['squiggle_start', 'squiggle_duration', 'po_width', 'po_height', 'po_x_width', 'po_y_min', 'po_y_max',
              'label_height']
toggle_inputs = ['toggle_x_axis', 'toggle_y_axis', 'toggle_annotations', 'toggle_smoothing']

app_data['wdg_dict'] = init_wdg_dict()
app_data['controls'] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg['wdg_width']))
p = figure(
    toolbar_location=None
)
p.outline_line_color = None
p.toolbar.logo = None
p.xaxis.visible = False
p.yaxis.visible = False
app_data['pore_plt'] = p

layout = row(
    app_data['controls'],
    app_data['pore_plt']
)

curdoc().title = "bulkvis"
curdoc().add_root(layout)
