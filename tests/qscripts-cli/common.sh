#
# Common stuff all the test scripts share.
#
# This script should be sourced in each run_test.sh.
#


TEST_ROOT=$(pwd)
test_name=$(basename ${TEST_ROOT})
# https://unix.stackexchange.com/questions/30091/fix-or-alternative-for-mktemp-in-os-x
TEST_TMP="$(mktemp -q -d -t ${test_name}.XXXXXX 2>/dev/null || mktemp -q -d)"

DIFFS="diffs.txt"
STDOUT="stdout.txt"
OK="\033[0;32mOK \033[0m"
FAIL="
\033[0;31mFAIL \033[0m Check the following files:
${TEST_TMP}/${DIFFS} 
${TEST_TMP}/${STDOUT}
"

# evaluates a string argument (diff command)
# if output code is 0: OK, if 1+: FAIL and log output
# example argument:
# diff <(egrep -v "version" test/probr_cl.fep) <(egrep -v "version" output/probr_cl.fep)
run_diff () {
  d=$(eval $@)
  if [[ $? -ne 0 ]]; then
    echo -e "${FAIL}"
    if [[ ${TRAVIS} == "true" ]]; then
        echo "$@"
        echo "${d}"
        # since we can't access the files on travis, print them out...
        echo "Outputs:"
        cat ${STDOUT}
    else
        echo "$@"   >> ${DIFFS}
        echo "${d}" >> ${DIFFS}
    fi
    if [[ ${PASSONFAIL} -eq 0 ]]; then
        exit 1
    fi
  else
    echo -e "${OK}"
  fi
}


