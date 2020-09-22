#!/bin/bash

set -e
set -eo pipefail

python -m pytest --cov=./ --cov-report=xml tests/
