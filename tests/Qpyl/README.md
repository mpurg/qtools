#### Running the tests
  
- Set the QBIN_DIR environment variable to point to the directory
of your Q6 binaries.  
- Make sure your are testing the correct qtools (source qtools_init.sh).
- Install pytest and run it in this directory.

```
# (in Bash)
export QBIN_DIR="/home/asd/apps/Q6/bin"
source ../../qtools_init.sh

pip install pytest 
pytest
```

