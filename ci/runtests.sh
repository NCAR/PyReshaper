#!/bin/bash

set -e
set -eo pipefail

python -m unittest discover --start-directory=source/test/ --pattern='*Tests.py'
