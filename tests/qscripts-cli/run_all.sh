#!/bin/bash

ROOTDIR=$(pwd)

for i in */
do
    cd $i;
    for j in */
    do
        cd $j
        echo "####################################################"
        echo "# Testing:    $i$j"
        echo "####################################################"
        ./run_test.sh
        cd $ROOTDIR/$i
    done
    cd $ROOTDIR
done

