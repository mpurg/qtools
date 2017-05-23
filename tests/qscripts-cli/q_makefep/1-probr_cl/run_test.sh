#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh


echo "Working in ${TEST_TMP}"
cp -r input ${TEST_TMP}
cd ${TEST_TMP}


q_makefep.py -s ./input/probr_cl_start.pdb \
             -m ./input/probr_cl.qmap \
             -l ./input/0-ff/*lib \
             -p ./input/0-ff/*prm \
             -f oplsaa \
             -o probr_cl.fep >> ${STDOUT}

for i in probr_cl.fep
do
    echo "Checking '${i}'"
    run_diff "diff <(sed 's/#.*//' ${i}) <(sed 's/#.*//' ${TEST_ROOT}/output/${i})"
done

cd ${TEST_ROOT}
