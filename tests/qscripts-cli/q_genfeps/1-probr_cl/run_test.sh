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
             relax \
             --repeats 1 \
             --frames 51 \
             --fromlambda 0.5 \
             --pdb ./input/0-topol/probr_cl_start.pdb \
             --prefix test_ >> ${STDOUT}

for i in test_000/*
do
    sed -i 's/#.*//' $i
done

echo "Checking outputs"
run_diff "diff test_000 ${TEST_ROOT}/output/test_000"

cd ${TEST_ROOT}
