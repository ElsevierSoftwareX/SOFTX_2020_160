#!/bin/bash -l

cd $2
nohup ./run_ts.py >> $1 2>&1 &

