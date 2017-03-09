#!/bin/bash 
QBIN=$(which qdyn5_r8)

OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

steps=( $(ls -1v *inp | sed 's/.inp//') )

for step in ${steps[@]}
do
 echo "Running step ${step}"
 if ${QBIN} ${step}.inp > ${step}.log
 then echo -e "$OK"
   cp ${step}.re  ${step}.re.rest
 else 
  echo -e "$FAILED"
  echo "Check output (${step}.log) for more info."
  exit 1
 fi
done

