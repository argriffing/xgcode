Description
-----------

This project serves a small dynamic web application without
requiring root privileges.
Multiple web apps are made available through a basic plugin architecture.
The main purpose of this project is to consolidate a bunch of random scripts.


Requirements
------------

Operation system requirements:

This project was developed using Ubuntu,
so it will probably work on Debian-based Linux distributions.
It probably will not work on Windows.

Major dependencies:

* A recent version of [python 2](http://www.python.org/) (2.6+).
* The [argparse](http://code.google.com/p/argparse/) module.
* The [cherrypy](http://www.cherrypy.org/) module (3.0+).
* The [epydoc](http://epydoc.sourceforge.net/) documentation program.
* The [cairo](http://www.cairographics.org/pycairo/) graphics library.
* The [numpy](http://numpy.scipy.org/) matrix and numerical library.
* The [scipy](http://www.scipy.org/) scientific library.

Minor dependencies
(each required by one or more individual scripts):

* The [matplotlib](http://matplotlib.sourceforge.net/) matlab-like library.


Web interface usage
-------------------

First create an empty directory which will be filled
with various temporary files created by the server.
Make this empty directory the current directory.

To tell the web app to start listening on 127.0.0.1:8080:

    python /path/to/run.py

For more options:

    python /path/to/run.py --help

You might want to run the program [in the background]
(http://apps.carleton.edu/curricular/cs/resources/source/bg_progs/).


Command line interface usage
----------------------------

Some of the individual scripts do useful things
when run by themselves on the command line.
To see which options are available for this run mode, try:

    python /path/to/20yymmdd.py --help

A script can be profiled using the command line modification:

    python -m cProfile /path/to/20yymmdd.py
