#!/bin/bash

set -e
set -eo pipefail

python -m pytest tests/
