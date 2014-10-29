==========
PyReshaper
==========

A package for converting NetCDF files from time-slice (history) format 
to time-series (single-variable) format.

:AUTHORS: John Dennis, Sheri Mickelson, Kevin Paul, Haiying Xu
:VERSION: 0.9.1
:COPYRIGHT: See the document entitled LICENSE.txt

Send questions and comments to Kevin Paul (kpaul@ucar.edu).


OVERVIEW
========

The PyReshaper package is a Python-based package for performing time-slice
to time-series convertion of NetCDF files, compliant with the CF 1.6 
Conventions.  The PyReshaper package is designed to run in parallel to
maximize performance, with the parallelism implemented over variables
(i.e., task parallelism).  This means that the maximum parallelism
achieveable for a given operation is one core/processor per variables in
the time-slice NetCDF files.


DEPENDENCIES
============

The PyReshaper directly depends upon the PyNIO and mpi4py packages.  Access
and manipulation of the NetCDF files is done through PyNIO, and the parallelism
is implimented directly with mpi4py.  Implicit dependencies exists, as PyNIO
has its own dependencies (netCDF, NCL, and numpy) as does mpi4py (numpy and 
MPI).

Currently the explicit dependencies are known to be:

* PyNIO (>=1.4.1)
* mpi4py (>=1.3)

This implies the dependencies:

* PyNIO depends upon numpy (>=1.4), NCL, and netCDF
* mpi4py depends on numpy (>=1.4) and MPI
    
Additionally, the entire package is designed to work with Python v2.6 and up
to (but not including) Python v3.0.
 
The version requirements have not been rigidly tested, so earlier versions
may actually work.  No version requirement is made during installation, though,
so problems might occur if an earlier versions of these packages have been
installed.


OBTAINING THE SOURCE CODE
=========================

Currently, the most up-to-date source code is available via svn from the site::

    https://subversion.ucar.edu/asap/pyReshaper/tags/v0.9.0

The source is available in read-only mode to everyone, but special permissions
can be given to those to make changes to the source.


BUILDING & INSTALLATION
=======================

Installation of the PyReshaper is very simple.  After checking out the source
from the above svn link, via::

    svn co https://subversion.ucar.edu/asap/pyReshaper/tags/v0.9.0 pyReshaper

change into the top-level source directory and run the Python distutils
setup.  On unix, this involves::

    $  cd pyReshaper
    $  python setup.py install [--prefix=/path/to/install/location]
    
The prefix is optional, as the default prefix is typically /usr/local on
linux machines.  However, you must have permissions to write to the prefix
location, so you may want to choose a prefix location where you have write
permissions.  Like most distutils installations, you can alternatively
install the pyReshaper with the --user option, which will automatically
select (and create if it does not exist) the $HOME/.local directory in which
to install.  To do this, type (on unix machines)::

    $  python setup.py install --user
    
This can be handy since the site-packages directory will be common for all
user installs, and therefore only needs to be added to the PYTHONPATH once.

To install API documentation for developer use, you must run doxygen with
the command (on unix machines)::

    $  doxygen Doxyfile

The resulting API documentation will be placed in the apidocs/ directory.


BEFORE USING THE PYRESHAPER PACKAGE
===================================

Before the PyReshaper package can be used, you must make sure that the 
site-packages directory containing the 'pyreshaper' source directory is in
your PYTHONPATH.  Depending on the PREFIX used during installation, this
path will be::

    $PREFIX/lib/python2.X/site-packages

where X will be 6 or 7 (or other) depending on the version of Python that you
are using to install the package.

To use the PyReshaper scripts (namely, 'slice2series'), you must add the
script binary directory to your PATH.  Depending on the PREFIX used during
installation, this path will be::

    $PREFIX/bin/
    
Once the script binary directory has been added to your PATH and the 
site-packages directory has been added to your PYTHONPATH, you may use the
PyReshaper package without issue.


INSTRUCTIONS & USE
==================

For instructions on how to use the PyReshaper, see the additional documents
found in the apidocs/ and docs/ directories.

If you are a developer wanting to use the PyReshaper API directly from your
own Python code, please read the 'BUILDING & INSTALLATION' section above
for instructions on how to build the API documentation.  Once built, you
will be able to open the 'apidocs/index.html' page in any browser.

The docs/ directory contains user manual describing how to use the binary 
scripts from the command-line as well as how to use the PyReshaper from 
within Python.  Both this README and the User Manual are written in 
reStructuredText (ReST), and can easily be converted to HTML or many other
formats with the help of a tool such as Docutils or Sphinx.

