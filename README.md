⇜ bulkvis ⇝
============

An app written in Python3 using [Bokeh](https://github.com/bokeh/bokeh/) to visualise raw squiggle data from Oxford Nanopore Technologies (ONT) bulkfiles. 

Table of contents
=================

 * [bulkvis](#-bulkvis-)
 * [Table of Contents](#table-of-contents)
 * [Installation](#installation)
 * [Configuration](#configuration)
 * [Launch the app](#launch-the-app)
 * [Tutorial](#tutorial)
 * [Support](#support)

Installation
============

```bash
# Make a python3 virtual environment
$ mkdir ~/envs
$ cd ~/envs
$ python3 -m venv ~/envs/bulkvis

# Activate virtual environment
$ source ~/envs/bulkvis/bin/activate

# Clone the repo to your installation/projects directory
$ git clone https://github.com/LooseLab/bulkvis.git

# Install dependencies via pip
$ pip install -r requirements.txt
```

Configuration
=============
bulkvis uses some user defined parameter to set default plot options. These can change depending on your preferences.
To set the parameters:

```INI
[plot_opts]
wdg_width = 300
plot_width = 970
plot_height = 800
y_min = 0
y_max = 2200
x_width = 15
```

Launch the app
==============

```bash
$ bokeh serve --show vis
```

Tutorial
========

:wrench: Coming soon 🤫

Support
=======

Open an issue please
