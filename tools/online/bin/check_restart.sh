#!/bin/bash -l

#alarm() { perl -e 'alarm shift; exec @ARGV' "$@"; }

dir=$1
log=$2
opt=$3
ex=$4

if [ ! -e $log ]; then
    touch $log
fi

cd $dir
echo $ex $log $dir $opt

if [ "$SITE_CLUSTER" = "ATLAS" ]; then
   /usr/bin/run-one $ex $log $dir $opt
else
   $CWB_ONLINE/bin/check_restart.py $log $ex $dir $opt
fi
