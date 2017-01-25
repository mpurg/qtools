#!/bin/bash

# make sure we are using this copy of the scripts
source ../../../../qtools_init.sh

# get the diffing function and some variable
source ../../common.sh

ROOTDIR=$(pwd)
mkdir -p $TESTDIR
cd $TESTDIR

q_makefep.py -s ../input/probr_cl_start.pdb \
             -m ../input/probr_cl.qmap \
             -l ../input/0-ff/*lib \
             -p ../input/0-ff/*prm \
             -f oplsaa \
             -o probr_cl.fep >> ${STDOUT}

for i in probr_cl.fep
do
    echo "Checking '$i'"
    run_diff "diff <(sed 's/#.*//' $i) <(sed 's/#.*//' ../output/$i)"
done

cd $ROOTDIR
