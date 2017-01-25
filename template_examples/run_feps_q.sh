#!/bin/bash 
CORES=4

QBIN=$(which qdyn5p)

OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

steps=( $(ls -1v *inp | sed 's/.inp//') )

for step in ${steps[@]}
do
 echo "Running step ${step}"
 if mpirun --report-bindings -n ${CORES} ${QBIN} ${step}.inp > ${step}.log
 then echo -e "$OK"
   cp ${step}.re  ${step}.re.rest
 else 
  echo -e "$FAILED"
  echo "Check output (${step}.log) for more info."
  exit 1
 fi
done

