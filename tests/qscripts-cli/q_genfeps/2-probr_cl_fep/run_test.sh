#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh


echo "Working in ${TEST_TMP}"
cp -r input ${TEST_TMP}
cd ${TEST_TMP}


q_genfeps.py ./input/genfeps.proc \
             ./input/1-relax/relax_003.inp \
             inp \
             --repeats 3 \
             --frames 11 \
             --fep ./input/0-topol/probr_cl.fep \
             --pdb ./input/0-topol/probr_cl_start.pdb \
             --prefix test_ >> ${STDOUT}

for i in test_*/*inp test_*/*fep
do
    sed -e 's/#.*//' -e 's/random_seed.*//' $i > ${i}.tmp
    mv ${i}.tmp $i
done

echo "Checking outputs"
run_diff "diff test_000 ${TEST_ROOT}/output/test_000"
run_diff "diff test_001 ${TEST_ROOT}/output/test_001"
run_diff "diff test_002 ${TEST_ROOT}/output/test_002"

cd ${TEST_ROOT}
