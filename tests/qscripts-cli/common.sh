#
# common stuff all the test scripts share
#

PASSONFAIL=0
if [[ $1 == "pass" ]]; then
    PASSONFAIL=1
fi

TESTDIR="tmp-run"
DIFFS="diffs.txt"
STDOUT="stdout.txt"
OK="\033[0;32mOK \033[0m"
FAIL="
\033[0;31mFAIL \033[0m Check the following files:
${TESTDIR}/${DIFFS} 
${TESTDIR}/${STDOUT}
"

# evaluates a string argument (diff command)
# if output code is 0: OK, if 1+: FAIL and log output
# example argument:
# diff <(egrep -v "version" test/probr_cl.fep) <(egrep -v "version" output/probr_cl.fep)
run_diff () {
  d=$(eval $@)
  if [[ $? -ne 0 ]]; then
    echo "$@"   >> ${DIFFS}
    echo "${d}" >> ${DIFFS}
    echo -e "${FAIL}"
    if [[ ${PASSONFAIL} -eq 0 ]]; then
        exit 1
    fi
  else
    echo -e "${OK}"
  fi
}


