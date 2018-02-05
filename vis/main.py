from collections import OrderedDict
import configparser
import math
import os
import re

import h5py
import numpy as np
import pandas as pd
from bokeh.plotting import curdoc, figure
from bokeh.layouts import row, widgetbox
from bokeh.models import TextInput, Button, Toggle, Div, Range1d, Label, Span, CheckboxGroup, Dropdown, PreText

config = configparser.ConfigParser()
config.read(os.path.dirname(os.path.realpath(__file__)) + '/config.ini')
cfg = config['plot_opts']


def init_wdg_dict():
    wdg_dict = OrderedDict()
    wdg_dict['data_input'] = TextInput(
        title='Data Source:',
        value="/Users/Alex/projects/bokehApp/bulkvis/data/NA12878_Run6_Sample2_29123.fast5",
        placeholder="Path to your bulkfile"
    )
    wdg_dict['data_button'] = Button(label="Update source file", button_type="primary")
    wdg_dict['data_button'].on_click(update_file)
    return wdg_dict


def update_file():
    """"""
    if app_data['bulkfile']:
        app_data['bulkfile'].flush()
        app_data['bulkfile'].close()

    file_src = app_data['wdg_dict']['data_input'].value
    # Clear old bulkfile data and build new data structures
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
            # get dataset length in seconds
            app_data['app_vars']['len_ds'] = math.ceil(len(signal_ds) / app_data['app_vars']['sf'])

    # add fastq and position inputs
    app_data['wdg_dict'] = init_wdg_dict()
    app_data['wdg_dict']['data_input'].value = file_src
    app_data['wdg_dict']['position'] = TextInput(
        title="Position:",
        value="",
        placeholder="Enter position or FASTQ ID",
        css_classes=[]
    )

    app_data['wdg_dict']['position'].on_change("value", parse_position)

    layout.children[0] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg['wdg_width']))


def open_bulkfile(path):
    """"""
    # !!! add in check to see if this is a ONT bulkfile
    # Open bulkfile in read-only mode
    file = h5py.File(path, "r")
    # Get sample frequency, how many data points are collected each second
    sf = int(file["UniqueGlobalKey"]["context_tags"].attrs["sample_frequency"].decode('utf8'))
    # make channel_list
    channel_list = np.arange(1, len(file["Raw"]) + 1, 1).tolist()
    return file, sf, channel_list


# noinspection PyUnboundLocalVariable
def parse_position(attr, old, new):
    if new[0] == "@":
        fq = new[1:]
        fq_list = fq.split(" ")
        for k, item in enumerate(fq_list):
            if k == 0:
                read_id = item
            if item.split("=")[0] == "ch":
                channel_num = item.split("=")[1]
                channel_str = "Channel_{num}".format(num=channel_num)
        # Get ch_str, start, end
        # If read_id and ch not set...
        # noinspection PyUnboundLocalVariable
        if read_id and channel_str:
            int_data_path = app_data['bulkfile']["IntermediateData"][channel_str]["Reads"]
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
            start_time = math.floor(df.iloc[0, :].read_start)
            end_time = math.ceil(df.iloc[-1, :].read_start)
        app_data['wdg_dict']['position'].value = "{ch}:{start}-{end}".format(ch=channel_num, start=start_time, end=end_time)
    elif re.match(r'^[0-9]{1,3}:[0-9]{1,9}-[0-9]{1,9}', new):
        coords = new.split(":")
        times = coords[1].split("-")
        channel_num = coords[0]
        channel_str = "Channel_{num}".format(num=channel_num)
        (start_time, end_time) = times[0], times[1]
    else:
        channel_str = None
        channel_num = None
        start_time = None
        end_time = None


    app_data['app_vars']['channel_str'] = channel_str
    app_data['app_vars']['channel_num'] = int(channel_num)
    app_data['app_vars']['start_time'] = int(start_time)
    app_data['app_vars']['end_time'] = int(end_time)
    update()


def update_data(bulkfile, app_vars):
    app_vars['duration'] = app_vars['end_time'] - app_vars['start_time']
    # get times and squiggles
    app_vars['start_squiggle'] = math.floor(app_vars['start_time'] * app_vars['sf'])
    app_vars['end_squiggle'] = math.floor(app_vars['end_time'] * app_vars['sf'])
    # get data in numpy arrays
    step = 1 / app_vars['sf']
    app_data['x_data'] = np.arange(app_vars['start_time'], app_vars['end_time'], step)
    app_data['y_data'] = bulkfile["Raw"][app_vars['channel_str']]["Signal"][()]
    app_vars['len_ds'] = len(app_data['y_data']) / app_vars['sf']
    app_data['y_data'] = app_data['y_data'][app_vars['start_squiggle']:app_vars['end_squiggle']]
    # get annotations
    path = bulkfile["IntermediateData"][app_vars['channel_str']]["Reads"]
    fields = ['read_id', 'read_start', 'modal_classification']
    app_data['int_label_df'], app_data['label_dt'] = get_annotations(path, fields)

    app_data['int_label_df'] = app_data['int_label_df'].drop_duplicates(subset="read_id", keep="first")
    app_data['int_label_df'].read_start = app_data['int_label_df'].read_start / app_vars['sf']
    app_data['int_label_df'].read_id = app_data['int_label_df'].read_id.str.decode('utf8')

    path = bulkfile["StateData"][app_vars['channel_str']]["States"]
    fields = ['analysis_raw_index', 'summary_state']
    app_data['state_label_df'], state_label_dtypes = get_annotations(path, fields)
    app_data['label_dt'].update(state_label_dtypes)


def get_annotations(path, fields):
    """
    Given a path to a hdf5 data set and a list of fields, return a
    :param path: path to hdf5 dataset eg: file[path][to][dataset]
    :param fields: list of fields in the dataset eg: ['list', 'of', 'fields']
    :return: pandas df of fields and dictionary of enumerated
    """
    data_labels = {}
    for field in fields:
        if h5py.check_dtype(enum=path.dtype[field]):
            dataset_dtype = h5py.check_dtype(enum=path.dtype[field])
            # data_dtype may lose some dataset dtypes there are duplicates of 'v'
            data_dtypes = {v: k for k, v in dataset_dtype.items()}
        data_labels[field] = path[field]

    labels_df = pd.DataFrame(data=data_labels)
    return labels_df, data_dtypes


def update():
    update_data(
        app_data['bulkfile'],
        app_data['app_vars']
    )
    if app_data['INIT']:
        build_widgets()
        layout.children[0] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg['wdg_width']))
        app_data['INIT'] = False
    app_data['wdg_dict']['duration'].text = "Duration: {d} seconds".format(d=app_data['app_vars']['duration'])
    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data,
        app_data['wdg_dict'],
        app_data['app_vars']
    )


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
    wdg['duration'] = PreText(text="Duration: {d} seconds".format(d=app_data['app_vars']['duration']))
    wdg['jump_next'] = Dropdown(label="Jump to next", button_type="primary", menu=jump_list)
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

    # wdg['label_filter'].on_change('active', update)
    wdg['jump_next'].on_click(next_update)
    # wdg['jump_prev'].on_click(prev_update)

    for name in toggle_inputs:
        wdg[name].on_click(toggle_button)
    for name in int_inputs:
        wdg[name].on_change('value', is_input_int)

    return wdg


def create_figure(x_data, y_data, app_data, wdg, app_vars):
    if wdg["toggle_smoothing"].active:
        w_range = app_vars['duration']
        divisor = math.e ** 2.5
        thin_factor = math.ceil(w_range / divisor)
    else:
        thin_factor = 1
    if thin_factor == 0:
        thin_factor = 1

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
            app_vars['start_time'], app_vars['start_time'] + int(wdg['po_x_width'].value)
        )
    p.xaxis.major_label_orientation = math.radians(45)

    if wdg['toggle_annotations'].active:
        slim_label_df = app_data['int_label_df'][
            (app_data['int_label_df']['read_start'] >= app_vars['start_time']) &
            (app_data['int_label_df']['read_start'] <= app_vars['end_time'])
        ]

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
                    text=str(app_data['label_dt'][label.modal_classification]),
                    level='glyph',
                    x_offset=0,
                    y_offset=0,
                    render_mode='canvas',
                    angle=-300
                )
                p.add_layout(labels)
        slim_label_df = app_data['state_label_df'][
            (app_data['state_label_df']['analysis_raw_index'] >= app_vars['start_time']) &
            (app_data['state_label_df']['analysis_raw_index'] <= app_vars['end_time'])
        ]

        for index, label in slim_label_df.iterrows():
            if app_data['label_mp'][label.summary_state] in wdg['label_filter'].active:
                event_line = Span(
                    location=label.analysis_raw_index,
                    dimension='height',
                    line_color='green',
                    line_dash='dashed',
                    line_width=1
                )
                p.add_layout(event_line)
                labels = Label(
                    x=label.analysis_raw_index,
                    y=int(wdg['label_height'].value),
                    text=str(app_data['label_dt'][label.summary_state]),
                    level='glyph',
                    x_offset=0,
                    y_offset=0,
                    render_mode='canvas',
                    angle=-300
                )
                p.add_layout(labels)
    return p


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
    update()


def toggle_button(state):
    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data,
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
    jump_start = app_data['int_label_df'][
        (app_data['int_label_df']['read_start'] > app_data['app_vars']['start_time'] + 1) &
        (app_data['int_label_df']['modal_classification'] == value)
    ]
    app_data['app_vars']['start_time'] = int(math.floor(jump_start['read_start'].iloc[0]))
    app_data['app_vars']['end_time'] = app_data['app_vars']['start_time'] + app_data['app_vars']['duration']
    app_data['wdg_dict']['position'].value = "{ch}:{start}-{end}".format(
        ch=app_data['app_vars']['channel_num'],
        start=app_data['app_vars']['start_time'],
        end=app_data['app_vars']['end_time']
    )
    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data,
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


app_data = {
    'file_src': None,  # bulkfile path (string)
    'bulkfile': None,  # bulkfile object
    'x_data': None,  # numpy ndarray time points
    'y_data': None,  # numpy ndarray signal data
    'int_label_df': None,  # pandas df of signal labels
    'state_label_df': None,  # pandas df of signal labels
    'label_dt': None,  # dict of signal enumeration
    'label_mp': None,  # dict matching labels to widget filter
    'app_vars': {  # dict of variables used in plots and widgets
        'len_ds': None,  # length of signal dataset in (seconds or events)
        'start_time': None,  # squiggle start time in seconds
        'end_time': None,  # squiggle end time in seconds
        'duration': None,  # squiggle duration in seconds
        'start_squiggle': None,  # squiggle start position (events)
        'end_squiggle': None,  # squiggle end position (events)
        'channel_str': None,  # 'Channel_NNN' (string)
        'channel_num': None,  # Channel number (int)
        'sf': None,  # sample frequency (int)
        'channel_list': None,  # list of all channels as int
        # 'fastq_str': None,              # fastq sequence id (string)
        # 'fastq_read_id': None           # fastq read id (string)
    },
    'wdg_dict': None,  # dictionary of widgets
    'controls': None,  # widgets added to widgetbox
    'pore_plt': None,  # the squiggle plot
    'INIT': True  # Initial plot with bulkfile (bool)
}

int_inputs = ['po_width', 'po_height', 'po_x_width', 'po_y_min', 'po_y_max', 'label_height']
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

curdoc().add_root(layout)
curdoc().title = "bulkvis"
