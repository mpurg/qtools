

OK="(\033[0;32m   OK   \033[0m)"
FAILED="(\033[0;31m FAILED \033[0m)"

steps=( $(ls -1v *inp | sed 's/.inp//') )

for step in ${steps[@]}
do
  echo "Running equilibration step ${step}"
  if qdyn5_r8 ${step}.inp > ${step}.log
  then
    echo -e "$OK"
  else 
    echo -e "$FAILED"
    echo "Check output (${step}.log) for more info."
    exit 1
  fi
done
