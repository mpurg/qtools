#!/bin/bash

source ../qtools_init.sh
echo -e "\nRunning tests of Qtools residing in \"${QTOOLS_HOME}\"\n"
sleep 2

if [[ -x "${QBIN_DIR}/Qdyn6" ]] && [[ -x "${QBIN_DIR}/Qfep6" ]] && [[ -x "${QBIN_DIR}/Qcalc6" ]]; then
    export QBIN_DIR="$( cd "${QBIN_DIR}" && pwd )" # convert to absolute
    echo "Q6 binaries found in ${QBIN_DIR}"
else
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
fi

echo
echo "--------------------------------------------------------------"
echo "|                     Running Qpyl tests                     |"
echo "--------------------------------------------------------------"

cd Qpyl
py.test --cov Qpyl || exit 1

echo
echo "--------------------------------------------------------------"
echo "|                     Running CLI tests                      |"
echo "--------------------------------------------------------------"

cd ../qscripts-cli
./run_all.sh || exit 1
