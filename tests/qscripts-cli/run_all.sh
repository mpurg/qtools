#!/bin/bash


export PASSONFAIL=0
if [[ $1 == "pass" ]]; then
    export PASSONFAIL=1
fi

export QTOOLS_TMPDIR=$(mktemp -d --tmpdir QTOOLS_TESTS_XXX)
TESTS_ROOT=$(pwd)

# make sure we are running this copy of the scripts
source ../../qtools_init.sh

echo "Setting up the CLI"
# backup existing config files
mv -f ${QTOOLS_HOME}/qscripts.cfg ${QTOOLS_HOME}/qscripts.cfg.bak
# run the config
qscripts_config.py --no_bin 2>&1

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
            mv -f ${QTOOLS_HOME}/qscripts.cfg.bak ${QTOOLS_HOME}/qscripts.cfg
            exit 1
        fi
    done
done

mv -f ${QTOOLS_HOME}/qscripts.cfg.bak ${QTOOLS_HOME}/qscripts.cfg
exit 0

