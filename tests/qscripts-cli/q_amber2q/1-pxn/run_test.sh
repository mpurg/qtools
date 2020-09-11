#!/bin/bash

# get the diffing function and some variable
source ../../common.sh

# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh

echo "Working in ${TEST_TMP}"
cp -r input ${TEST_TMP}
cd ${TEST_TMP}

q_amber2q.py -m ./input/pxn.mol2 \
             -i ./input/pxn.prepi \
             -p ./input/gaff2_at16.dat \
             -f ./input/pxn_gaff2.frcmod \
             --ignore_errors \
             -o pxn >> ${STDOUT}

echo ${TESTS_ROOT}
for i in pxn.lib pxn.prm pxn.prm.chk
do
    echo "Checking '${i}'"
    run_diff "diff <(sed 's/#.*//' ${i}) <(sed 's/#.*//' ${TEST_ROOT}/output/${i})"
done

cd ${TEST_ROOT}
