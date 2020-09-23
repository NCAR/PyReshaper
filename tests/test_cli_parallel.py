"""
Copyright 2020, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

from glob import glob
from os import remove
from os.path import dirname, join, realpath
from subprocess import Popen

import pytest

from .checks import check_outfile
from .data import config

THISDIR = dirname(realpath(__file__))
INFILES = sorted(glob(join(THISDIR, 'data', 'input*.nc')))
METAFILE = join(THISDIR, 'data', 'metafile.nc')


def cmd_opts(**kwargs):
    args = kwargs.pop('args', [])
    opts = []
    for opt in kwargs:
        opt_str = '--{}'.format(opt)
        opt_vals = kwargs[opt] if isinstance(kwargs[opt], list) else [kwargs[opt]]
        for opt_val in opt_vals:
            if isinstance(opt_val, bool):
                if opt_val:
                    opts.append(opt_str)
            else:
                opts.extend([opt_str, str(opt_val)])
    opts.extend(args)
    return opts


def check(**make_args):
    kwargs = {
        'infiles': make_args.pop('args', []),
        'prefix': make_args.pop('output_prefix', 'tseries.'),
        'suffix': make_args.pop('output_suffix', '.nc'),
        'metadata': make_args.pop('metadata', []),
        'once': make_args.pop('once', False),
    }
    kwargs.update(make_args)
    return check_outfile(**kwargs)


@pytest.mark.parametrize('n', [0, 1, 2, 4])
@pytest.mark.parametrize(
    'make_args',
    [
        {
            'args': INFILES,
        },
        {
            'args': INFILES,
            'compression_level': 0,
            'output_prefix': 'out.',
            'metadata': [v for v in config.tvmvars] + ['time'] + [v for v in config.chvars],
        },
    ],
)
@pytest.mark.parametrize(
    'run_args',
    [
        {
            'args': ['input.s2s'],
        },
    ],
)
def test_cli_mpi(n, make_args, run_args):
    pre_cmds = ['coverage', 'run', '-p', '-m']

    spec_cmds = pre_cmds + ['pyreshaper.cli.s2smake', '-o', 'input.s2s'] + cmd_opts(**make_args)
    p_spec = Popen(spec_cmds)
    p_spec.communicate()
    assert p_spec.returncode == 0

    mpirun = ['mpirun', '-n', str(n)] if n > 0 else []
    run_cmds = mpirun + pre_cmds + ['pyreshaper.cli.s2srun'] + cmd_opts(**run_args)
    p_run = Popen(run_cmds)
    p_run.communicate()
    assert p_run.returncode == 0

    for tsvar in config.tsvars:
        check(tsvar=tsvar, **make_args)

    remove('input.s2s')
    for fname in glob('*.nc'):
        remove(fname)
