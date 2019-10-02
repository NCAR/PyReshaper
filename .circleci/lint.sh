#!/bin/bash

set -e
set -eo pipefail

echo "Code Styling with flake8 and isort"

source activate ${ENV_NAME}

echo "[flake8]"
flake8

echo "[isort]"
isort --recursive -w 100 --check-only source
