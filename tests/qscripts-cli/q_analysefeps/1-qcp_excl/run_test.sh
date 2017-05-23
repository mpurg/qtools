#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh


echo "Working in ${TEST_TMP}"
cp -r input ${TEST_TMP}
cd ${TEST_TMP}

q_analysefeps.py ./input ./input3 \
                 --qfep_out qfep_123.log \
                 --out qaf.out \
                 --plots_out qaf.pd.json \
                 --subcalcs --subcalc_dir random.dir_name \
                 --lra_l 0.9 0.4   >> ${STDOUT}

for i in qaf.out qaf.pd.json random.dir_name/qaf*json
do
    echo "Checking '${i}'"
    run_diff "diff <(sed 's/#.*//' ${i}) <(sed 's/#.*//' ${TEST_ROOT}/output/${i})"
done

cd ${TEST_ROOT}
