#!/bin/bash

source ../qtools_init.sh
echo -e "\nRunning tests of Qtools residing in \"${QTOOLS_HOME}\"\n\n"

cd Qpyl
py.test --cov Qpyl

cd ../qscripts-cli
./run_all.sh
