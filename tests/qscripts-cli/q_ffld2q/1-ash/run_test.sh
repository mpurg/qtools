#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh



ROOTDIR=$(pwd)
mkdir -p $TESTDIR
cd $TESTDIR
rm ace_ash_nma.lib ace_ash_nma.prm ace_ash_nma.prm.chk 2>/dev/null

q_ffld2q.py ../input/ace_ash_nma.ffld11 ../input/ace_ash_nma.pdb -o ace_ash_nma >> ${STDOUT}

for i in ace_ash_nma.lib ace_ash_nma.prm ace_ash_nma.prm.chk
do
    echo "Checking '$i'"
    run_diff "diff <(sed 's/#.*//' $i) <(sed 's/#.*//' ../output/$i)"
done

cd $ROOTDIR
