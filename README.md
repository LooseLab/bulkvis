⇜ bulkvis ⇝
============

An app written in Python3 using [Bokeh][1] to visualise raw squiggle data from Oxford Nanopore Technologies (ONT) bulkfiles. 

Quickstart
==========
We require python 3.6 or newer.

Using `conda` with this environment setup:
```yaml
name: bulkvis
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.7
  - pip
  - pip:
    - git+https://github.com/LooseLab/bulkvis.git@2.0
```

or with another python source:

```bash
# Make a python3 virtual environment
python3 -m venv bulkvis

# Activate virtual environment
source bulkvis/bin/activate

# Clone the repo to your installation/projects directory
pip install git+https://github.com/LooseLab/bulkvis.git@2.0

# Start bokeh server
bulkvis serve <BULK_FILE_DIRECTORY> --show
```

Other install requires:
===

To open some bulk FAST5 files [`vbz compression plugins`][2] are required. 
These are written and maintained by Oxford Nanopore Technologies.


 [1]: https://github.com/bokeh/bokeh/
 [2]: https://github.com/nanoporetech/vbz_compression
