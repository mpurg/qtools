#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh


echo "Working in ${TEST_TMP}"
cp -r input ${TEST_TMP}
cd ${TEST_TMP}


q_ffld2q.py ./input/ace_ash_nma.ffld14 ./input/ace_ash_nma.pdb -o ace_ash_nma >> ${STDOUT}

for i in ace_ash_nma.lib ace_ash_nma.prm ace_ash_nma.prm.chk
do
    echo "Checking '$i'"
    run_diff "diff <(sed 's/#.*//' $i) <(sed 's/#.*//' ${TEST_ROOT}/output/$i)"
done

cd ${TEST_ROOT}
