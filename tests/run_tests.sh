#!/bin/bash

source ../qtools_init.sh
echo -e "\nRunning tests of Qtools residing in \"${QTOOLS_HOME}\"\n"
sleep 2

echo "Downloading Q6 binaries"
mkdir -p qbin && cd qbin
export QBIN_DIR=$(pwd)
echo $QBIN_DIR
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    wget -q -O- https://github.com/qusers/Q6/releases/download/v6.0/Q6Linux.tar.gz | tar xv
elif [[ "$OSTYPE" == "darwin"* ]]; then
    wget -q -O- https://github.com/qusers/Q6/releases/download/v6.0/Q6Mac.tar.gz | tar xv
else
    echo "Precompiled binaries for ${OSTYPE} are not available. Some tests will fail."
    sleep 3
fi
cd ..

cd Qpyl
py.test --cov Qpyl

cd ../qscripts-cli
./run_all.sh
