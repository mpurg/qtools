#!/bin/bash

# find the directory where this script is located in
export QTOOLS_HOME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# set the path for the CLI scripts
export PATH=${QTOOLS_HOME}/qscripts-cli:${PATH}

# set the python path for packages (Qpif)
export PYTHONPATH=${QTOOLS_HOME}/packages:${PYTHONPATH}

