#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh


echo "Working in ${TEST_TMP}"
cp -r ../../shared_data/probr_cl/* ${TEST_TMP}
cd ${TEST_TMP}


q_automapper.py 12 -5 0 0 \
               --dirs testrep1 testrep2 \
               --skip 1 \
               --min 1 \
               --bins 50 \
               --temp 298 \
               --step 10 \
               --threshold 0.001 \
               --iter 10 \
               --nosingle \
               --qfep_exec ${QBIN_DIR}/Qfep6 \
               --out qautomapper.log > qautomapper1.out

q_automapper.py 12 -5 0 0 \
               --dirs testrep1 testrep2 \
               --skip 1 \
               --min 1 \
               --bins 50 \
               --temp 298 \
               --step 10 \
               --threshold 10 \
               --iter 10 \
               --nosingle \
               --qfep_exec ${QBIN_DIR}/Qfep6 \
               --out asd > qautomapper2.out

q_automapper.py 12 -5 0 0 \
               --dirs testrep1 testrep2 \
               --skip 1 \
               --min 1 \
               --bins 50 \
               --temp 298 \
               --step 1 \
               --threshold 0.001 \
               --iter 10 \
               --nosingle \
               --qfep_exec ${QBIN_DIR}/Qfep6 \
               --out asd > qautomapper3.out

q_automapper.py 12 -5 -100 0 \
               --dirs testrep1 testrep2 \
               --skip 1 \
               --min 1 \
               --bins 50 \
               --temp 298 \
               --step 1 \
               --threshold 0.001 \
               --iter 5 \
               --nosingle \
               --qfep_exec ${QBIN_DIR}/Qfep6 \
               --out asd > qautomapper4.out

q_automapper.py 12 15 0 0 \
               --dirs testrep1 testrep2 \
               --skip 1 \
               --min 1 \
               --bins 50 \
               --temp 298 \
               --step 10 \
               --threshold 0.001 \
               --iter 10 \
               --nosingle \
               --qfep_exec ${QBIN_DIR}/Qfep6 \
               --out asd > qautomapper5.out

for i in *log *out
do
    echo "Checking '${i}'"
    run_diff "diff <(sed 's/#.*//' ${i}) <(sed 's/#.*//' ${TEST_ROOT}/output/${i})"
done

cd ${TEST_ROOT}
