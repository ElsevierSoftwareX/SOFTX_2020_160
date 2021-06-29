#!/bin/bash -l

cd $2
nohup ./run.py >> $1 2>&1 &

