#!/bin/bash

set -e
set -eo pipefail

source activate ${ENV_NAME}

echo "[Serial Tests]"
python -m unittest discover --start-directory=source/test/ --pattern='*Tests.py'
