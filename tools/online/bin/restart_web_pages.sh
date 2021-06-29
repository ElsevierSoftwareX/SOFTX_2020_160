#!/bin/bash -l

cd $2
nohup ${CWB_ONLINE}/bin/web_pages.sh $2 $3 > $1 2>&1 &
