⇜ bulkvis ⇝
============

An app written in Python3 using [Bokeh](https://github.com/bokeh/bokeh/) to visualise 
raw squiggle data from Oxford Nanopore Technologies (ONT) bulkfiles. 

Quickstart
==========
```bash
# Make a python3 virtual environment
$ python3 -m venv bulkvis

# Activate virtual environment
$ source bulkvis/bin/activate

# Clone the repo to your installation/projects directory
$ pip install git+https://github.com/LooseLab/bulkvis.git@2.0

# Start bokeh server
$ bulkvis serve <BULK_FILE_DIRECTORY> --show
```

Other install requires:
===

To open some bulk FAST5 files [`vbz compression plugins`](https://github.com/nanoporetech/vbz_compression) 
are required. These are written and maintained by Oxford Nanopore Technologies.
