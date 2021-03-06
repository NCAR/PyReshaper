#!/usr/bin/env python
"""
This script provides the command-line interface (CLI) to the PyReshaper

This script is designed to run a specfile (i.e., a Pickled Specifier object).
The specfile itself should be constructed from a hand-written Python script,
or from the makes2sspec tool that accompanies this script.


Copyright 2020 University Corporation for Atmospheric Research

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import optparse
import pickle

from pyreshaper import reshaper


def cli(argv=None):
    desc = """This tool is designed to run a PyReshaper Specifier as read from a pickled Specifier object (specfile)"""

    parser = optparse.OptionParser(prog='s2srun', description=desc)
    parser.add_option(
        '-1',
        '--once',
        default=False,
        action='store_true',
        dest='once',
        help=('Whether to write a "once" file with all ' 'metadata. [Default: False]'),
    )
    parser.add_option(
        '-c',
        '--read_chunk',
        default=None,
        action='append',
        dest='rchunks',
        help=(
            'Read chunk size for a named dimension.  This should '
            'be given as a comma-separated pair (e.g., NAME,SIZE) '
            'indicating the name of the dimension to chunk over '
            'and the chunk size.  Multiple chunk options can be '
            'given on the command line, each one enabling chunking '
            'over a new dimension.  [Default: None]'
        ),
    )
    parser.add_option(
        '-l',
        '--limit',
        default=0,
        type='int',
        help=(
            'The limit on the number of time-series files per '
            'processor to write.  Useful when debugging.  A '
            'limit of 0 means write all output files.'
            '[Default: 0]'
        ),
    )
    parser.add_option(
        '-m',
        '--write_mode',
        default='w',
        type='str',
        help=(
            'Determine the write mode to use for writing '
            "output files.  Can be 'w' for normal operation, "
            "'s' to skip generation of time-series files if "
            "the files already exist, 'o' to overwrite "
            "existing time-series files, or 'a' to append "
            "to existing time-series files [Default: 'w']"
        ),
    )
    parser.add_option(
        '-s',
        '--serial',
        default=False,
        action='store_true',
        dest='serial',
        help=('Whether to run in serial (True) or parallel ' '(False). [Default: False]'),
    )
    parser.add_option(
        '-v',
        '--verbosity',
        default=1,
        type='int',
        help=(
            'Verbosity level for level of output.  A value '
            'of 0 means no output, and a value greater than '
            '0 means more output detail. [Default: 1]'
        ),
    )
    parser.add_option(
        '-w',
        '--write_chunk',
        default=None,
        action='append',
        dest='wchunks',
        help=(
            'Write chunk size for a named dimension, used when '
            'writing time-series variables only (default chunking '
            'is used for metadata variables).  This should '
            'be given as a comma-separated pair (e.g., NAME,SIZE) '
            'indicating the name of the dimension to chunk over '
            'and the chunk size.  Multiple chunk options can be '
            'given on the command line, each one enabling chunking '
            'over a new dimension.  [Default: None]'
        ),
    )

    opts, args = parser.parse_args(argv)

    if len(args) == 0:
        raise ValueError('Must supply a specfile as input')
    else:
        specfile = args[0]

    if opts.rchunks is not None:
        opts.rchunks = dict((c.split(',')[0], int(c.split(',')[1])) for c in opts.rchunks)

    if opts.wchunks is not None:
        opts.wchunks = dict((c.split(',')[0], int(c.split(',')[1])) for c in opts.wchunks)

    return opts, specfile


def main(argv=None):
    opts, specfile = cli(argv)

    # Try importing the file
    try:
        spec = pickle.load(open(specfile, 'rb'))
    except:
        err_msg = "Specifier File '{}' could not be opened and read".format(specfile)
        raise RuntimeError(err_msg)

    # Create the PyReshaper object
    reshpr = reshaper.create_reshaper(
        spec, serial=opts.serial, verbosity=opts.verbosity, wmode=opts.write_mode, once=opts.once
    )

    # Run the conversion (slice-to-series) process
    reshpr.convert(output_limit=opts.limit, rchunks=opts.rchunks, wchunks=opts.wchunks)

    # Print timing diagnostics
    reshpr.print_diagnostics()


if __name__ == '__main__':
    main()
