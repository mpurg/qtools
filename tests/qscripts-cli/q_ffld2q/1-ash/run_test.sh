#!/bin/bash

# make sure we are using this copy of the scripts
source ../../../../qscripts_init.sh

# get the diffing function and some variable
source ../../common.sh

ROOTDIR=$(pwd)
mkdir -p $TESTDIR
cd $TESTDIR

q_ffld2q.py ../input/ace_ash_nma.ffld11 ../input/ace_ash_nma.pdb -o ace_ash_nma >> ${STDOUT}

for i in ace_ash_nma.lib ace_ash_nma.prm ace_ash_nma.prm.chk
do
    echo "Checking '$i'"
    run_diff "diff <(sed 's/#.*//' $i) <(sed 's/#.*//' ../output/$i)"
done

cd $ROOTDIR
