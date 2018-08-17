import configparser
from dateutil import parser
import math
from pathlib import Path
import re
from collections import OrderedDict

import h5py
import numpy as np
import pandas as pd
from bokeh.layouts import row, widgetbox
from bokeh.models import TextInput, Toggle, Div, Range1d, Label, Span, Title, LabelSet, RadioButtonGroup
from bokeh.models import CheckboxGroup, Dropdown, PreText, Select, Button, ColumnDataSource
from bokeh.plotting import curdoc, figure
from utils.stitch import export_read_file

config = configparser.ConfigParser()
config.read(str(Path(Path(__file__).resolve().parent / 'config.ini')))
cfg_po = config['plot_opts']
cfg_dr = config['data']
cfg_lo = config['labels']
output_backend = {'canvas', 'svg', 'webgl'}


"""

PUT IN A CHECK THAT REQUIRED CFG PARAMS ARE SET!!!!

if cfg_dr[out] == '':
    disable read-file write function <- remove button

"""


def init_wdg_dict():
    """
    Initialise the widget dictionary, adds in the initial bulkfile selector
    Returns
    -------

    """
    wdg_dict = OrderedDict()
    wdg_dict['file_list'] = Select(title="Select bulk-file:", options=app_data['app_vars']['files'])
    wdg_dict['file_list'].on_change('value', update_file)
    return wdg_dict


def update_file(attr, old, new):
    """"""
    if app_data['bulkfile']:
        app_data['bulkfile'].flush()
        app_data['bulkfile'].close()

    if new == "":
        app_data['wdg_dict'] = init_wdg_dict()
        f = figure(toolbar_location=None)
        f.outline_line_color = None
        f.toolbar.logo = None
        f.xaxis.visible = False
        f.yaxis.visible = False
        layout.children[0] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg_po['wdg_width']))
        layout.children[1] = f
        return

    file_src = app_data['wdg_dict']['file_list'].value
    file_wdg = app_data['wdg_dict']['file_list']
    file_list = app_data['app_vars']['files']
    map_file_list = app_data['app_vars']['map_files']
    # Clear old bulkfile data and build new data structures
    app_data.clear()
    app_data['app_vars'] = {}
    app_data['wdg_dict'] = OrderedDict()
    app_data['label_dt'] = OrderedDict()
    app_data['file_src'] = Path(Path(cfg_dr['dir']) / file_src)
    app_data['INIT'] = True
    app_data['app_vars']['files'] = file_list
    app_data['app_vars']['map_files'] = map_file_list

    (app_data['bulkfile'],
     app_data['app_vars']['sf'],
     app_data['app_vars']['attributes']) = open_bulkfile(app_data['file_src'])

    raw_path = app_data['bulkfile']["Raw"]
    for i, member in enumerate(raw_path):
        if i == 0:
            signal_ds = raw_path[member]["Signal"][()]
            # get dataset length in seconds
            # app_data['app_vars']['len_ds'] = math.ceil(len(signal_ds) / app_data['app_vars']['sf'])
            app_data['app_vars']['len_ds'] = len(signal_ds) / app_data['app_vars']['sf']

    # add fastq and position inputs
    app_data['wdg_dict'] = init_wdg_dict()
    app_data['wdg_dict']['file_list'] = file_wdg
    # app_data['wdg_dict']['maps_list'] = Select(title="Select mapping file:", options=map_file_list)
    app_data['wdg_dict']['position_label'] = Div(text='Position', css_classes=['position-dropdown', 'help-text'])
    app_data['wdg_dict']['position_text'] = Div(
        text="""Enter a position in your bulkfile as <code>channel:start_time-end_time</code> or a
                <code>complete FASTQ header</code>.
                """,
        css_classes=['position-drop']
    )
    app_data['wdg_dict']['position'] = TextInput(
        value="",
        placeholder="e.g 391:120-150 or complete FASTQ header",
        css_classes=['position-label']
    )
    read_bmf(app_data['app_vars']['Run ID'])
    app_data['wdg_dict']['position'].on_change("value", parse_position)

    layout.children[0] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg_po['wdg_width']))


def read_bmf(run_id):
    run_id = run_id + '.bmf'
    try:
        app_data['bmf'] = pd.read_csv(Path(Path(cfg_dr['map']) / run_id), sep='\t')
        # filter mappings to just this run
        app_data['bmf'] = app_data['bmf'][app_data['bmf']['run_id'] == app_data['app_vars']['Run ID']]
    except FileNotFoundError:
        pass
    except Exception as e:
        print(e)
    return


def open_bulkfile(path):
    # !!! add in check to see if this is a ONT bulkfile
    # Open bulkfile in read-only mode
    open_file = h5py.File(path, "r")
    # Get sample frequency, how many data points are collected each second
    sf = int(open_file["UniqueGlobalKey"]["context_tags"].attrs["sample_frequency"].decode('utf8'))
    attributes = OrderedDict([
        ('tracking_id', [('Experiment', 'sample_id'), ('Flowcell ID', 'flow_cell_id'), ('MinKNOW version', 'version'),
                         ('Protocols version', 'protocols_version'), ('MinION ID', 'device_id'),
                         ('Hostname', 'hostname'), ('Run ID', 'run_id'), ('ASIC ID', 'asic_id'),
                         ('Experiment start', 'exp_start_time')]),
        ('context_tags', [('Sequencing kit', 'sequencing_kit'), ('Flowcell type', 'flowcell_type')])])

    for k, v in attributes.items():
        for attribute in v:
            try:
                app_data['app_vars'][attribute[0]] = open_file['UniqueGlobalKey'][k].attrs[attribute[1]].decode('utf8')
                if attribute[1] == 'exp_start_time':
                    app_data['app_vars'][attribute[0]] = parser.parse(
                        app_data['app_vars'][attribute[0]]).strftime('%d-%b-%Y %H:%M:%S')
            except KeyError:
                app_data['app_vars'][attribute[0]] = 'N/A'
    return open_file, sf, attributes


# noinspection PyUnboundLocalVariable
def parse_position(attr, old, new):
    if re.match(r'^(\@[a-f0-9\-]{36})([a-z0-9=\s]{1,})ch=[0-9]{1,4}', new):
        # https://regex101.com/r/9VvgNM/4
        # Match UUID / read_id as fastq str
        #   ^(\@[a-f0-9\-]{36})
        # Match lowercase a-z, 0-9, '=' and whitespace
        #   ([a-z0-9=\s]{1,})
        # Match 'ch=' and up to 4 numbers
        #   ch=[0-9]{1,4}
        # if new[0] == "@":
        input_error(app_data['wdg_dict']['position'], 'remove')
        fq = new[1:]
        fq_list = fq.split(" ")
        # split out read_id and channel
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
            if len(df) > 2:
                start_time = math.floor(df.iloc[0, :].read_start)
                end_time = math.ceil(df.iloc[-1, :].read_start)
            else:
                input_error(app_data['wdg_dict']['position'], 'add')
                return
        else:
            input_error(app_data['wdg_dict']['position'], 'add')
            return
        app_data['wdg_dict']['position'].value = "{ch}:{start}-{end}".format(
            ch=channel_num,
            start=start_time,
            end=end_time
        )
    elif re.match(r'^([0-9]{1,4}:[0-9]{1,9}-[0-9]{1,9})\Z', new):
        # https://regex101.com/r/zkN1j2/2
        input_error(app_data['wdg_dict']['position'], 'remove')
        coords = new.split(":")
        times = coords[1].split("-")
        channel_num = coords[0]
        channel_str = "Channel_{num}".format(num=channel_num)
        (start_time, end_time) = int(times[0]), int(times[1])
        if end_time - start_time <= 0:
            input_error(app_data['wdg_dict']['position'], 'add')
            return
    else:
        input_error(app_data['wdg_dict']['position'], 'add')
        return

    if int(end_time) > app_data['app_vars']['len_ds']:
        end_time = app_data['app_vars']['len_ds']
    app_data['app_vars']['channel_str'] = channel_str
    app_data['app_vars']['channel_num'] = int(channel_num)
    app_data['app_vars']['start_time'] = int(start_time)
    app_data['app_vars']['end_time'] = int(end_time)

    app_data['wdg_dict']['position'].value = "{ch}:{start}-{end}".format(
        ch=app_data['app_vars']['channel_num'],
        start=app_data['app_vars']['start_time'],
        end=app_data['app_vars']['end_time']
    )

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
    app_data['label_df'], app_data['label_dt'] = get_annotations(path, fields, 'modal_classification')
    app_data['label_df'] = app_data['label_df'].drop_duplicates(subset=['read_id', 'modal_classification'], keep="first")
    app_data['label_df'].read_start = app_data['label_df'].read_start / app_vars['sf']
    app_data['label_df'].read_id = app_data['label_df'].read_id.str.decode('utf8')

    path = bulkfile["StateData"][app_vars['channel_str']]["States"]
    fields = ['acquisition_raw_index', 'summary_state']
    state_label_df, state_label_dtypes = get_annotations(path, fields, 'summary_state')
    state_label_df.acquisition_raw_index = state_label_df.acquisition_raw_index / app_vars['sf']
    state_label_df = state_label_df.rename(
        columns={'acquisition_raw_index': 'read_start', 'summary_state': 'modal_classification'}
    )
    app_data['label_df'] = app_data['label_df'].append(state_label_df, ignore_index=True)
    app_data['label_df'].sort_values(by='read_start', ascending=True, inplace=True)
    app_data['label_dt'].update(state_label_dtypes)


def get_annotations(path, fields, enum_field):
    data_labels = {}
    for field in fields:
        data_labels[field] = path[field]
    data_dtypes = {}
    if h5py.check_dtype(enum=path.dtype[enum_field]):
        dataset_dtype = h5py.check_dtype(enum=path.dtype[enum_field])
        # data_dtype may lose some dataset dtypes there are duplicates of 'v'
        data_dtypes = {v: k for k, v in dataset_dtype.items()}
    labels_df = pd.DataFrame(data=data_labels)
    return labels_df, data_dtypes


def build_widgets():
    """"""
    check_labels = []
    jump_list = []
    check_active = []
    app_data['label_mp'] = {}
    for k, v in enumerate(app_data['label_dt'].items()):
        app_data['label_mp'][v[0]] = k
        check_labels.append(v[1])
        if v[1] in cfg_lo:
            if cfg_lo[v[1]] == 'True':
                check_active.append(k)
                jump_list.append((v[1], str(v[0])))
        else:
            print("label {v} is in your bulk-file but not defined in config.ini".format(v=v[1]))
            check_active.append(k)

    if len(check_active) == len(check_labels):
        filter_toggle_active = 0
    elif len(check_active) == 0:
        filter_toggle_active = 1
    else:
        filter_toggle_active = None

    wdg = app_data['wdg_dict']
    wdg['duration'] = PreText(text="Duration: {d} seconds".format(d=app_data['app_vars']['duration']), css_classes=['duration_pre'])
    wdg['navigation_label'] = Div(text='Navigation:', css_classes=['navigation-dropdown', 'help-text'])
    wdg['navigation_text'] = Div(
        text="""Use the <code><b>Jump to ...</b></code> buttons to find the next or previous event type.
                """,
        css_classes=['navigation-drop']
    )
    wdg['jump_next'] = Dropdown(label="Jump to next", button_type="primary", menu=jump_list, css_classes=['jump-block'])
    wdg['jump_prev'] = Dropdown(label="Jump to previous", button_type="primary", menu=jump_list)

    wdg['export_label'] = Div(text='Export data:', css_classes=['export-dropdown', 'help-text'])
    wdg['export_text'] = Div(
        text="""Export data, as a read file, from the current position. These are written to the output directory 
                specified in your config file.
                """,
        css_classes=['export-drop']
    )
    wdg['save_read_file'] = Button(
        label="Save read file",
        button_type="success",
        css_classes=[]
    )
    wdg['bulkfile_info'] = Div(text='Bulkfile info', css_classes=['bulkfile-dropdown', 'caret-down'])
    wdg['bulkfile_help'] = Div(text='Bulkfile info help:', css_classes=['bulkfile-help-dropdown', 'help-text', 'bulkfile-drop'])
    wdg['bulkfile_help_text'] = Div(
        text="""This contains basic information about the experiment that is recorded in the bulk-fast5-file.
                """,
        css_classes=['bulkfile-help-drop']
    )
    wdg['bulkfile_text'] = Div(text="", css_classes=['bulkfile-drop'])
    for k, v in app_data['app_vars']['attributes'].items():
        for entry in v:
            wdg['bulkfile_text'].text += '<b>{f}:</b> <br><code>{val}</code><br>'.format(
                f=entry[0],
                val=app_data['app_vars'][entry[0]]
            )
    wdg['label_options'] = Div(text='Select annotations', css_classes=['filter-dropdown', 'caret-down'])
    wdg['filter_help'] = Div(text='filter help:', css_classes=['filter-help-dropdown', 'help-text', 'filter-drop'])
    wdg['filter_help_text'] = Div(
        text="""Select which bulkfile annotations should be rendered on the chart. 'Display annotations' will turn all 
                annotations on or off.
                """,
        css_classes=['filter-help-drop']
    )
    wdg['toggle_annotations'] = Toggle(
        label="Display annotations",
        button_type="danger",
        css_classes=['toggle_button_g_r', 'filter-drop'],
        active=True
    )
    wdg['toggle_mappings'] = Toggle(
        label="Display mappings",
        button_type="danger",
        css_classes=['toggle_button_g_r', 'filter-drop'],
        active=True
    )
    wdg['filter_toggle_group'] = RadioButtonGroup(labels=["Select all", "Select none"], active=filter_toggle_active, css_classes=['filter-drop'])
    wdg['label_filter'] = CheckboxGroup(labels=check_labels, active=check_active, css_classes=['filter-drop'])
    
    wdg['plot_options'] = Div(text='Plot adjustments', css_classes=['adjust-dropdown', 'caret-down'])
    wdg['adjust_help'] = Div(text='adjust help:', css_classes=['adjust-help-dropdown', 'help-text', 'adjust-drop'])
    wdg['adjust_help_text'] = Div(
        text="""Adjust chart parameters, such as width, height and where annotations are rendered. These are set in the
                config.ini, where the default values can be edited.
                """,
        css_classes=['adjust-help-drop']
    )
    wdg['po_width'] = TextInput(title='Plot Width (px)', value=cfg_po['plot_width'], css_classes=['adjust-drop'])
    wdg['po_height'] = TextInput(title='Plot Height (px)', value=cfg_po['plot_height'], css_classes=['adjust-drop'])
    wdg['label_height'] = TextInput(
        title="Annotation height (y-axis)",
        value=cfg_po['label_height'],
        css_classes=['adjust-drop']
    )
    wdg['po_y_max'] = TextInput(title="y max", value=cfg_po['y_max'], css_classes=['adjust-drop', 'toggle_y_target'])
    wdg['po_y_min'] = TextInput(title="y min", value=cfg_po['y_min'], css_classes=['adjust-drop', 'toggle_y_target'])
    wdg['toggle_y_axis'] = Toggle(
        label="Fixed Y-axis",
        button_type="danger",
        css_classes=['toggle_button_g_r', 'adjust-drop', 'toggle_y_axis'],
        active=False
    )
    wdg['toggle_smoothing'] = Toggle(
        label="Smoothing",
        button_type="danger",
        css_classes=['toggle_button_g_r', 'adjust-drop'],
        active=True
    )

    wdg['label_filter'].on_change('active', update_checkboxes)
    wdg['filter_toggle_group'].on_change('active', update_toggle)
    wdg['jump_next'].on_click(next_update)
    wdg['jump_prev'].on_click(prev_update)
    wdg['save_read_file'].on_click(export_data)

    for name in toggle_inputs:
        wdg[name].on_click(toggle_button)
    for name in int_inputs:
        wdg[name].on_change('value', is_input_int)
    return wdg


def create_figure(x_data, y_data, wdg, app_vars):


    def vline(x_coords, y_upper, y_lower):
        # Return a dataset that can plot vertical lines
        x_values = np.vstack((x_coords, x_coords)).T
        y_upper_list = np.full((1, len(x_values)), y_upper)
        y_lower_list = np.full((1, len(x_values)), y_lower)
        y_values = np.vstack((y_lower_list, y_upper_list)).T
        return x_values.tolist(), y_values.tolist()


    def hlines(y_coords, x_lower, x_upper):
        """

        Parameters
        ----------
        y_coords: (int, float) height to plot lines at
        x_lower: (int, float) lower x coord
        x_upper: (int, float) upper x coord

        Returns
        -------

        """
        x_values = np.vstack((x_lower, x_upper)).T
        y_values_list = np.full((1, len(x_values)), y_coords)
        y_values = np.vstack((y_values_list, y_values_list)).T
        return x_values.tolist(), y_values.tolist()
    
    
    if wdg["toggle_smoothing"].active:
        w_range = app_vars['duration']
        divisor = math.e ** 2.5
        thin_factor = math.ceil(w_range / divisor)
    else:
        thin_factor = 1
    if thin_factor == 0:
        thin_factor = 1

    greater_delete_index = np.argwhere(y_data > int(cfg_po['upper_cut_off']))
    x_data = np.delete(x_data, greater_delete_index)
    y_data = np.delete(y_data, greater_delete_index)

    lesser_delete_index = np.argwhere(y_data < int(cfg_po['lower_cut_off']))
    x_data = np.delete(x_data, lesser_delete_index)
    y_data = np.delete(y_data, lesser_delete_index)

    x_data = x_data[::thin_factor]
    y_data = y_data[::thin_factor]

    data = {
        'x': x_data,
        'y': y_data,
    }

    source = ColumnDataSource(data=data)

    p = figure(
        plot_height=int(wdg['po_height'].value),
        plot_width=int(wdg['po_width'].value),
        toolbar_location="right",
        tools=['xpan', 'xbox_zoom', 'undo', 'reset', 'save'],
        active_drag="xpan",
    )
    if cfg_po['output_backend'] not in output_backend:
        p.output_backend = 'canvas'
    else:
        p.output_backend = cfg_po['output_backend']
    # Add step/% points plotted: Step: {sp} ({pt:.3f}) -> sp=thin_factor, pt=1/thin_factor
    p.add_layout(Title(
        text="Channel: {ch} Start: {st} End: {ed} Sample rate: {sf}".format(
            ch=app_vars['channel_num'],
            st=app_vars['start_time'],
            ed=app_vars['end_time'],
            sf=app_vars['sf']
        )),
        'above'
    )
    p.add_layout(Title(
        text="Bulk-file: {s}".format(s=app_data['wdg_dict']["file_list"].value)),
        'above'
    )

    p.toolbar.logo = None
    p.yaxis.axis_label = "Raw signal"
    p.yaxis.major_label_orientation = "horizontal"
    p.xaxis.axis_label = "Time (seconds)"
    p.line(source=source, x='x', y='y', line_width=1)
    p.xaxis.major_label_orientation = math.radians(45)
    p.x_range.range_padding = 0.01

    # set padding manually
    y_min = np.amin(data['y'])
    y_max = np.amax(data['y'])
    pad = (y_max - y_min)*0.1 / 2
    p.y_range = Range1d(y_min - pad, y_max + pad)
    try:
        app_data['bmf']
    except NameError:
        bmf_set = False
    except KeyError:
        bmf_set = False
    else:
        bmf_set = True
    if bmf_set and wdg['toggle_mappings'].active:
        # set padding manually
        # lower_pad = (y_max - y_min) * 0.1 / 2
        # upper_pad = (y_max - y_min) / 2
        # p.y_range = Range1d(y_min - lower_pad, y_max + upper_pad)
        # set mapping track midpoints
        # upper_mapping = upper_pad / 4 * 3 + y_max
        # lower_mapping = upper_pad / 4 + y_max
        lower_mapping = int(wdg['label_height'].value) + 750
        # Select only this channel
        slim_bmf = app_data['bmf'][app_data['bmf']['channel'] == app_vars['channel_num']]
        # Select the current viewed range
        slim_bmf = slim_bmf[
            (
                    (slim_bmf['start_time'] > app_vars['start_time']) &
                    (slim_bmf['end_time'] < app_vars['end_time'])
            ) |
            (
                    (slim_bmf['start_time'] < app_vars['start_time']) &
                    (slim_bmf['end_time'] < app_vars['end_time']) &
                    (slim_bmf['end_time'] > app_vars['start_time'])
            ) |
            (
                    (slim_bmf['start_time'] > app_vars['start_time']) &
                    (slim_bmf['end_time'] > app_vars['end_time']) &
                    (slim_bmf['start_time'] < app_vars['end_time'])
            )
        ]
        slim_bmf['start_time'] = slim_bmf['start_time'].where(slim_bmf['start_time'] > app_vars['start_time'], app_vars['start_time'])
        slim_bmf['end_time'] = slim_bmf['end_time'].where(slim_bmf['end_time'] < app_vars['end_time'], app_vars['end_time'])

        slim_bmf['height'] = lower_mapping
        slim_bmf['offset'] = np.ones(len(slim_bmf)) * 5 + slim_bmf.groupby(['start_time', 'end_time']).cumcount() * 15
        # Convert slim_bmf to ColDataSrc
        mapping_source = ColumnDataSource(data=slim_bmf.to_dict(orient='list'))
        # Add labels to LabelSet
        mapping_labels = LabelSet(x='start_time', y='height', text='label', level='glyph', x_offset=5,
                                  y_offset='offset', source=mapping_source, render_mode='canvas')
        p.add_layout(mapping_labels)
        # Add some colour here
        # Forward Vertical lines => blue
        p_x, p_y = vline(np.concatenate([slim_bmf['start_time'].where(slim_bmf['strand'] == '+').dropna().values,
                                         slim_bmf['end_time'].where(slim_bmf['strand'] == '+').dropna().values]),
                         lower_mapping + 20,
                         lower_mapping - 20
                         )
        p.multi_line(p_x, p_y, line_dash='solid', color='blue', line_width=1)
        # Reverse Vertical lines => red
        p_x, p_y = vline(np.concatenate([slim_bmf['start_time'].where(slim_bmf['strand'] == '-').dropna().values,
                                         slim_bmf['end_time'].where(slim_bmf['strand'] == '-').dropna().values]),
                         lower_mapping + 20,
                         lower_mapping - 20
                         )
        p.multi_line(p_x, p_y, line_dash='solid', color='red', line_width=1)
        # Horizontal lines
        p_x, p_y = hlines(lower_mapping,
                          slim_bmf['start_time'].where(slim_bmf['strand'] == '+').dropna(),
                          slim_bmf['end_time'].where(slim_bmf['strand'] == '+').dropna()
        )
        p.multi_line(p_x, p_y, line_dash='solid', color='blue', line_width=1)
        # Horizontal lines
        p_x, p_y = hlines(lower_mapping,
                          slim_bmf['start_time'].where(slim_bmf['strand'] == '-').dropna(),
                          slim_bmf['end_time'].where(slim_bmf['strand'] == '-').dropna()
        )
        p.multi_line(p_x, p_y, line_dash='solid', color='red', line_width=1)
    
    if wdg['toggle_y_axis'].active:
        p.y_range = Range1d(int(wdg['po_y_min'].value), int(wdg['po_y_max'].value))
    if wdg['toggle_annotations'].active:
        # Map modal_classifications onto df
        app_data['label_df']['mc_active_map'] = app_data['label_df']['modal_classification'].map(app_data['label_mp'])
        app_data['label_df']['mc_label_map'] = app_data['label_df']['modal_classification'].map(app_data['label_dt'])
        # Here labels are thinned out
        slim_label_df = app_data['label_df'][
            (app_data['label_df']['read_start'] >= app_vars['start_time']) &
            (app_data['label_df']['read_start'] <= app_vars['end_time'])
            ]
        # Use pd.isin to remove unwanted annotations from the slimmed df
        slim_label_df = slim_label_df[slim_label_df['mc_active_map'].isin(wdg['label_filter'].active) == True]
        # get coordinates and vstack them to produce [[x, x], [x, x]...]
        line_x_values = np.vstack((slim_label_df['read_start'].values, slim_label_df['read_start'].values)).T
        tmp_list = np.full((1, len(line_x_values)), -10000)
        line_y_values = np.vstack((tmp_list, tmp_list * -1)).T
        # Add all vertical lines as multi_line
        p.multi_line(line_x_values.tolist(), line_y_values.tolist(), line_dash='dashed', color='green', line_width=1)
        # combine series to form label
        slim_label_df['label'] = slim_label_df['mc_label_map'] + " - " + slim_label_df['read_id'].astype('str')
        # Create ColumnDataSource combining labels and coordinates
        label_source = ColumnDataSource(
            data=dict(
                x=slim_label_df['read_start'].values,
                y=np.full((len(slim_label_df), 1), int(wdg['label_height'].value)),
                t=slim_label_df['label'].values
            )
        )
        # Add all labels as a label set
        labels = LabelSet(x='x', y='y', text='t', level='glyph',
                          x_offset=0, y_offset=0, source=label_source,
                          render_mode='canvas', angle=-270, angle_units='deg')
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


def update():
    update_data(
        app_data['bulkfile'],
        app_data['app_vars']
    )
    if app_data['INIT']:
        build_widgets()
        layout.children[0] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg_po['wdg_width']))
        app_data['INIT'] = False
    app_data['wdg_dict']['duration'].text = "Duration: {d} seconds".format(d=app_data['app_vars']['duration'])
    app_data['wdg_dict']['toggle_smoothing'].active = True
    layout.children[1] = create_figure(
        app_data['x_data'],
        app_data['y_data'],
        app_data['wdg_dict'],
        app_data['app_vars']
    )


def update_other(attr, old, new):
    update()


def update_toggle(attr, old, new):
    if new == 0:
        app_data['wdg_dict']['label_filter'].active = list(np.arange(0, len(app_data['wdg_dict']['label_filter'].labels), 1))
    elif new == 1:
        app_data['wdg_dict']['label_filter'].active = []
    update()


def update_checkboxes(attr, old, new):
    if len(new) != len(app_data['wdg_dict']['label_filter'].labels) and len(new) != 0:
        app_data['wdg_dict']['filter_toggle_group'].active = None
    update()


def next_update(value):
    if value != 'reset':
        value = int(value)
        jump_start = app_data['label_df'][
            (app_data['label_df']['read_start'] > app_data['app_vars']['start_time'] + 1) &
            (app_data['label_df']['modal_classification'] == value)
            ]
        try:
            app_data['app_vars']['start_time'] = int(math.floor(jump_start['read_start'].iloc[0]))
        except IndexError:
            app_data['wdg_dict']['duration'].text += "\n{ev} event not found".format(ev=app_data['label_dt'][value])
            return
        except Exception as e:
            print(type(e))
            print(e)
        app_data['app_vars']['end_time'] = app_data['app_vars']['start_time'] + app_data['app_vars']['duration']
        app_data['wdg_dict']['position'].value = "{ch}:{start}-{end}".format(
            ch=app_data['app_vars']['channel_num'],
            start=app_data['app_vars']['start_time'],
            end=app_data['app_vars']['end_time']
        )
        layout.children[1] = create_figure(
            app_data['x_data'],
            app_data['y_data'],
            app_data['wdg_dict'],
            app_data['app_vars']
        )
        app_data['wdg_dict']['jump_next'].value = "reset"
    else:
        return


def prev_update(value):
    if value != 'reset':
        value = int(value)
        jump_start = app_data['label_df'][
            (app_data['label_df']['read_start'] < app_data['app_vars']['start_time']) &
            (app_data['label_df']['modal_classification'] == value)
            ]
        try:
            app_data['app_vars']['start_time'] = int(math.floor(jump_start['read_start'].iloc[-1]))
        except IndexError:
            app_data['wdg_dict']['duration'].text += "\n{ev} event not found".format(ev=app_data['label_dt'][value])
            return
        except Exception as e:
            print(type(e))
            print(e)
        app_data['app_vars']['end_time'] = app_data['app_vars']['start_time'] + app_data['app_vars']['duration']
        app_data['wdg_dict']['position'].value = "{ch}:{start}-{end}".format(
            ch=app_data['app_vars']['channel_num'],
            start=app_data['app_vars']['start_time'],
            end=app_data['app_vars']['end_time']
        )
        layout.children[1] = create_figure(
            app_data['x_data'],
            app_data['y_data'],
            app_data['wdg_dict'],
            app_data['app_vars']
        )
        app_data['wdg_dict']['jump_prev'].value = "reset"
    else:
        return


def export_data():
    try:
        start_val = math.floor(app_data['app_vars']['start'] * app_data['app_vars']['sf'])
        end_val = math.ceil(app_data['app_vars']['end'] * app_data['app_vars']['sf'])
    except KeyError:
        start_val = app_data['app_vars']['start_squiggle']
        end_val = app_data['app_vars']['end_squiggle']
    if export_read_file(
        app_data['app_vars']['channel_num'],
        start_val,
        end_val,
        app_data['bulkfile'],
        cfg_dr['out']
    ) == 0:
        app_data['wdg_dict']['duration'].text += "\nread file created"
    else:
        app_data['wdg_dict']['duration'].text += "\nError: read file not created"


app_data = {
    'file_src': None,  # bulkfile path (string)
    'bulkfile': None,  # bulkfile object
    'bmf': None,  # bmf dataframe
    'x_data': None,  # numpy ndarray time points
    'y_data': None,  # numpy ndarray signal data
    'label_df': None,  # pandas df of signal labels
    'label_dt': None,  # dict of signal enumeration
    'label_mp': None,  # dict matching labels to widget filter
    'app_vars': {  # dict of variables used in plots and widgets
        'len_ds': None,  # length of signal dataset
        'start_time': None,  # squiggle start time in seconds
        'end_time': None,  # squiggle end time in seconds
        'duration': None,  # squiggle duration in seconds
        'start_squiggle': None,  # squiggle start position (samples)
        'end_squiggle': None,  # squiggle end position (samples)
        'channel_str': None,  # 'Channel_NNN' (string)
        'channel_num': None,  # Channel number (int)
        'sf': None,  # sample frequency (int)
        'attributes': None  # OrderedDict of bulkfile attr info
    },
    'wdg_dict': None,  # dictionary of widgets
    'controls': None,  # widgets added to widgetbox
    'pore_plt': None,  # the squiggle plot
    'INIT': True  # Initial plot with bulkfile (bool)
}

int_inputs = ['po_width', 'po_height', 'po_y_min', 'po_y_max', 'label_height']
toggle_inputs = ['toggle_y_axis', 'toggle_annotations', 'toggle_mappings', 'toggle_smoothing']

app_data['app_vars']['files'] = []
p = Path(cfg_dr['dir'])
app_data['app_vars']['files'] = [(x.name, x.name) for x in p.iterdir() if x.suffix == '.fast5']
m = Path(cfg_dr['map'])
app_data['app_vars']['map_files'] = [(x.name, x.name) for x in m.iterdir() if x.suffix == '.bmf']
app_data['app_vars']['map_files'].insert(0, ("", "--"))
# check files are useable by h5py
for index, file in enumerate(app_data['app_vars']['files']):
    file = file[0]
    try:
        bulk_file = h5py.File(Path(Path(cfg_dr['dir']) / file), 'r')
    except OSError:
        app_data['app_vars']['files'][index] = None
        continue
    try:
        try_path = bulk_file["Raw"]
    except KeyError:
        app_data['app_vars']['files'][index] = None
        continue
    for i, channel in enumerate(try_path):
        if i == 0:
            try:
                try_path[channel]["Signal"][0]
            except KeyError:
                app_data['app_vars']['files'][index] = None
        break
    bulk_file.flush()
    bulk_file.close()
app_data['app_vars']['files'] = list(filter((None).__ne__, app_data['app_vars']['files']))
app_data['app_vars']['files'].insert(0, ("", "--"))

app_data['wdg_dict'] = init_wdg_dict()
app_data['controls'] = widgetbox(list(app_data['wdg_dict'].values()), width=int(cfg_po['wdg_width']))

f = figure(toolbar_location=None)
f.line(x=[0], y=[0])
f.outline_line_color = None
f.toolbar.logo = None
f.xaxis.visible = False
f.yaxis.visible = False
f.xgrid.visible = False
f.ygrid.visible = False
app_data['pore_plt'] = f

layout = row(
    app_data['controls'],
    app_data['pore_plt']
)

curdoc().add_root(layout)
curdoc().title = "bulkvis"
