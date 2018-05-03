⇜ bulkvis ⇝
============

An app written in Python3 using [Bokeh](https://github.com/bokeh/bokeh/) to visualise raw squiggle data from Oxford Nanopore Technologies (ONT) bulkfiles. 
See the documentation at [https://bulkvis.readthedocs.io](https://bulkvis.readthedocs.io)

Quickstart
==========
```bash
# Make a python3 virtual environment
$ mkdir ~/envs
$ cd ~/envs
$ python3 -m venv ~/envs/bulkvis

# Activate virtual environment
$ source ~/envs/bulkvis/bin/activate

# Clone the repo to your installation/projects directory
$ git clone https://github.com/LooseLab/bulkvis.git

# Enter the bulkvis folder
$ cd bulkvis

# Install dependencies via pip
$ pip install -r requirements.txt

# Set config with set_config.py
$ python utils/set_config.py -b <<bulkfile>> -i /path/to/bulkfile/directory -e /path/to/read/file/export/directory

# Move to bulkvis' parent folder
$ cd ..

# Start bokeh server
$ bokeh serve --show bulkvis
```

Configuration
=============
bulkvis uses some user defined parameter to set default plot options. 
These can change depending on your preferences. For a comprehensive overview
of the config file see [config.md](config.md)
