#!/bin/bash

set -e
set -eo pipefail

source activate ${ENV_NAME}
pytest --junitxml=test-reports/junit.xml --cov=./ --verbose
echo "[Upload coverage]"
codecov
