#!/bin/bash

# make sure we are using this copy of the scripts
source ../../../../qtools_init.sh

# get the diffing function and some variable
source ../../common.sh

ROOTDIR=$(pwd)
mkdir -p $TESTDIR
cd $TESTDIR
rm qaf.out qaf.pd.json qaf.pd2.json qaf.pd3.json 2>/dev/null

q_analysefeps.py ../input ../input3 \
                 --qfep_out qfep_123.log \
                 --out qaf.out \
                 --plots_out qaf.pd.json \
                 --po_qcp qaf.pd2.json \
                 --po_qcpm qaf.pd3.json \
                 --lra_l 0.9 0.4   >> ${STDOUT}

for i in qaf.out qaf.pd.json qaf.pd2.json qaf.pd3.json;
do
    echo "Checking '$i'"
    run_diff "diff <(sed 's/#.*//' $i) <(sed 's/#.*//' ../output/$i)"
done

cd $ROOTDIR
