#!/bin/bash

# make sure we are using this copy of the scripts
source ../../../../qtools_init.sh

# get the diffing function and some variable
source ../../common.sh

ROOTDIR=$(pwd)
mkdir -p $TESTDIR
cd $TESTDIR

rm test_relax/* 2>/dev/null
rmdir test_relax 2>/dev/null

q_genrelax.py ../input/genrelax.proc \
              --top ../input/0-topol/probr_cl.top \
              --pdb ../input/0-topol/probr_cl_start.pdb \
              --fep ../input/0-topol/probr_cl.fep \
              --rs ../input/run_relax_q.sh \
              --outdir test_relax >> ${STDOUT}

for i in test_relax/*
do
    sed -i 's/#.*//' $i
done

echo "Checking outputs"
run_diff "diff test_relax $ROOTDIR/output/test_relax"

cd $ROOTDIR
