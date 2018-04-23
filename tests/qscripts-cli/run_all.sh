#!/bin/bash


export PASSONFAIL=0
if [[ $1 == "pass" ]]; then
    export PASSONFAIL=1
fi

export QTOOLS_TMPDIR=$(mktemp -d --tmpdir QTOOLS_TESTS_XXX)
TESTS_ROOT=$(pwd)

# make sure we are running this copy of the scripts
source ../../qtools_init.sh
# run the config script
qscripts_config.py 2>&1 >/dev/null

for i in q_*/
do
    cd ${TESTS_ROOT}/${i}
    for j in */
    do
        cd ${TESTS_ROOT}/${i}/${j}
        echo "####################################################"
        echo "|  Testing:    ${i}${j}"
        echo "####################################################"
        ./run_test.sh
        if [[ $? != 0 ]]; then
            exit 1
        fi
    done
done


