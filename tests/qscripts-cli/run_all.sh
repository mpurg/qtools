#!/bin/bash

export PASSONFAIL=0
if [[ $1 == "pass" ]]; then
    export PASSONFAIL=1
fi

ROOTDIR=$(pwd)

# make sure we are running this copy of the scripts
source ../../qtools_init.sh
# run the config script
qscripts_config.py 2>&1 >/dev/null

for i in */
do
    cd $i;
    for j in */
    do
        cd $j
        echo "####################################################"
        echo "# Testing:    $i$j"
        echo "####################################################"
        ./run_test.sh
        if [[ $? != 0 ]]; then
            exit 1
        fi
        cd $ROOTDIR/$i
    done
    cd $ROOTDIR
done

