############
Installation
############

We recommend running bulkvis from within a python virtual environment so that there are no conflicts in dependencies.

Installing pip
==============

pip is most likely already installed, to find out run::

    pip --version

If pip is not installed, use the official
`get-pip.py <https://pip.pypa.io/en/stable/installing/#installing-with-get-pip-py>`_ script.

Create and activate a virtual environment
=========================================

For linux and MacOS::

    python3 -m venv bulkvis-env
    source bulkvis-env/bin/activate

For Windows::

    python3 -m venv bulkvis
    bulkvis\Scripts\activate

If the virtual environment is successfully activated the prefix ``(bulkvis)`` will be present.

Running ``deactivate`` will deactivate and exit the virtual environment

Clone bulkvis
=============

bulkvis can be retrieved by cloning the git repository::

    git clone https://github.com/LooseLab/bulkvis.git

or by navigating to `bulkvis <https://github.com/LooseLab/bulkvis.git>`_ and downloading an zip of the repository,
this will then need to be unzipped.

Installing dependencies
=======================

Once the repository is cloned or downloaded bulkvis' dependencies will need to be installed. This **must** be run from
within the virtual environment to prevent conflicts. Run::

    pip install -r bulkvis/requirements.txt

This will fetch and install all the required packages.

Creating config.ini
===================

bulkvis uses a configuration file, config.ini, to provide global variables that are required to

Starting bulkvis
================

To start bulkvis::

    bokeh serve --show bulkvis

