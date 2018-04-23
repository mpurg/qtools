#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh


echo "Working in ${TEST_TMP}"
cp -r ../../shared_data/probr_cl/* ${TEST_TMP}
cd ${TEST_TMP}


q_mapper.py 0 0 \
            --dirs testrep1 testrep2 testrep3 \
            --skip 1 \
            --min 1 \
            --bins 20 \
            --temp 298 \
            --qfep_exec ${QBIN_DIR}/Qfep6 \
            --out qmapper1.out >> ${STDOUT}

q_mapper.py 100 10 \
            --dirs testrep1 testrep2 testrep3 \
            --skip 100 \
            --min 1 \
            --bins 20 \
            --temp 298 \
            --qfep_exec ${QBIN_DIR}/Qfep6 \
            --out qmapper2.out >> ${STDOUT}

q_mapper.py 100 10 \
            --skip 1 \
            --min 1 \
            --bins 20 \
            --temp 298 \
            --qfep_exec ${QBIN_DIR}/Qfep6 \
            --out qmapper3.out >> ${STDOUT}

for i in qmapper*out
do
    echo "Checking '${i}'"
    run_diff "diff <(sed 's/#.*//' ${i}) <(sed 's/#.*//' ${TEST_ROOT}/output/${i})"
done

cd ${TEST_ROOT}
