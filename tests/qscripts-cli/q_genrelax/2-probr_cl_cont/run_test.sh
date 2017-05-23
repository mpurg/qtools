#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh


echo "Working in ${TEST_TMP}"
cp -r input ${TEST_TMP}
cd ${TEST_TMP}


q_genrelax.py ./input/genrelax.proc \
              --cont ./input/1-relax/relax_003.inp \
              --rest cont_inp \
              --pdb ./input/0-topol/probr_cl_start.pdb \
              --rs ./input/run_relax_q.sh \
              --outdir test_relax >> ${STDOUT}

for i in test_relax/*
do
    sed -i 's/#.*//' ${i}
done

echo "Checking outputs"
run_diff "diff test_relax ${TEST_ROOT}/output/test_relax"

cd ${TEST_ROOT}
