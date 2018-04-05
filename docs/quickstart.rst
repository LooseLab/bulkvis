##########
Quickstart
##########

This page will provide a quick overview of the bulkvis features and how to use them.

Start bulkvis
=============
From the directory containing bulkvis, run::

    bokeh serve --show bulkvis

This will start the bulkvis app and open the page in your default web browser. The page should look like this:

.. figure:: _static/01_initial.png
    :class: figure
    :alt: Screenshot of bulkvis showing a blank screen with a drop-down box in the top left corner

    Screenshot of bulkvis on initial load

Selecting a bulk-fast5-file
===========================
Use the drop-down to select a bulk-fast5-file from your supplied directory. Once the file is loaded the position can be
entered.

.. figure:: _static/02_position.png
    :class: figure
    :alt: Screenshot of bulkvis showing, in the top left corner, a drop-down box with a file selected and a text box labeled position

Selecting a position
====================
The position can be selected by either supplying coordinates or a fastq read header. The contents of this input is submitted
by clicking away from the text box or by pressing return/enter. If bulkvis cannot parse the input the text box will turn
red until valid input is detected.

After a position is entered bulkvis will completely load and the chart will be visible.

Using coordinates
-----------------
Coordinates refer to the channel, start time and end time. This is given in the format ``channel:start-end``. For
example to navigate to channel 42 and see the squiggle from 30 seconds to 90 seconds::

    42:30-90

Using a fastq read header
-------------------------
Alternatively, the position can be given as a fastq read header that is from the run associated with this bulk-fast5-file.
This can be copied and pasted into the text-box e.g::

    @b45a4b09-6f22-40f6-afd9-aa7fca8e89f3 runid=f9291b45b0c66faa77755e51738d193fcfafffc7 read=234 ch=391 start_time=2018-01-18T21:59:40Z


After entering valid input the chart and other elements will load:

.. figure:: _static/03_plot.png
    :class: figure
    :alt: Screenshot of bulkvis showing a 'squiggle' plot of raw nanopore signal and a left-hand sidebar containing information about the plot



Navigating
==========
The bulk-fast5-file can be navigated by jumping to the next or previous event and by using the xpan (|xpan_icon|) to drag
the plot along the x-axis or zoom (|zoom_icon|) to take a closer look at a section of the plot.

The jump to action is available for any event type that is listed as ``True`` in config.ini and is available even when the
event is not currently being displayed.


Bulkfile information
====================
The bulkfile information panel showcases information that is present in the bulk-fast5-file that is not necessarily
displayed in MinKNOW.

.. figure:: _static/04_sidebar.png
    :class: figure
    :alt: Screenshot of the sidebar from bulkvis, showing the file selections drop-down, position input, jump-to buttons, export button, information panel, and two hidden sections ('Select annotations' and 'Plot adjustments')

Screenshot of sidebar

Annotations
===========
How to use annotations ...

Plot adjustments
================
What plot adjustments can be made

Exporting read files
====================
How to export read files and where they are written

Exporting images
================
How to export images of plots, this is browser dependent...

.. |zoom_icon| image:: /_static/icons/zoom.png
    :height: 11pt
.. |xpan_icon| image:: /_static/icons/xpan.png
    :height: 11pt