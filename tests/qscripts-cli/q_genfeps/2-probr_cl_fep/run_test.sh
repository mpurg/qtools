#!/bin/bash

# get the diffing function and some variable
source ../../common.sh


# make sure we are running this copy of the scripts
source ../../../../qtools_init.sh



ROOTDIR=$(pwd)
mkdir -p $TESTDIR
cd $TESTDIR

rm test_*/* 2>/dev/null
rmdir test_* 2>/dev/null

q_genfeps.py ../input/genfeps.proc \
             ../input/1-relax/relax_003.inp \
             inp \
             --repeats 3 \
             --frames 11 \
             --fep ../input/0-topol/probr_cl.fep \
             --pdb ../input/0-topol/probr_cl_start.pdb \
             --prefix test_ >> ${STDOUT}

for i in test_*/*
do
    sed -i -e 's/#.*//' -e 's/random_seed.*//' $i
done

echo "Checking outputs"
run_diff "diff test_000 $ROOTDIR/output/test_000"
run_diff "diff test_001 $ROOTDIR/output/test_001"
run_diff "diff test_002 $ROOTDIR/output/test_002"

cd $ROOTDIR
