#!/bin/bash

set -e
set -eo pipefail

coverage run -p -m pytest tests/
coverage combine
coverage xml
