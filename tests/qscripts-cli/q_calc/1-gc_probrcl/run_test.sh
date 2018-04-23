#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh


echo "Working in ${TEST_TMP}"
cp -r ../../shared_data/dfpase_dfp/* ${TEST_TMP}
cd ${TEST_TMP}

q_calc.py gc dfpase_dfp_start.pdb 15 20 --iscale 1 \
                                        --dirs testrep1 testrep2 testrep3 \
                                        --lra_l 0.86 0.46 \
                                        --plots_out qgc1.pd.json \
                                        --out qc1.log \
                                        --pdbgc gc1.pdb \
                                        --qcalc_exec ${QBIN_DIR}/Qcalc6 \
                                        --writeout > ${STDOUT}

#q_calc.py gc dfpase_dfp_start.pdb 15 20 --iscale 4 \
#                                        --dirs testrep1 testrep2 testrep3 \
#                                        --lra_l 0.86 0.46 \
#                                        --plots_out qgc2.pd.json \
#                                        --out qc2.log \
#                                        --pdbgc gc2.pdb \
#                                        --qcalc_exec ${QBIN_DIR}/Qcalc6 \
#                                        --writeout > ${STDOUT}
#
#echo "4841 4842 4843 4844" > qmask.list
#q_calc.py gc dfpase_dfp_start.pdb 15 20 --iscale 1 \
#                                        --dirs testrep1 testrep2 testrep3 \
#                                        --lra_l 0.86 0.46 \
#                                        --qmask qmask.list \
#                                        --plots_out qgc3.pd.json \
#                                        --out qc3.log \
#                                        --pdbgc gc3.pdb \
#                                        --qcalc_exec ${QBIN_DIR}/Qcalc6 \
#                                        --writeout > ${STDOUT}
#

for i in *log *json *pdb
do
    echo "Checking '${i}'"
    run_diff "diff <(sed 's/#.*//' ${i}) <(sed 's/#.*//' ${TEST_ROOT}/output/${i})"
done

cd ${TEST_ROOT}
