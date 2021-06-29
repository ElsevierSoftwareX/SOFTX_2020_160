#!/bin/bash

if [[ -z $1 ]]; then
  echo ""
  echo "cwb_pegasus_submit.sh :  missing dag file"
  echo ""
  exit 1
else
  CWB_DAGFILE=$1                              # DAG FILE
  echo "input  dag file : " \'$CWB_DAGFILE\'
  CWB_DAXFILE="${CWB_DAGFILE/.dag/.dax}"      # DAX FILE
  echo "output dax file : " \'$CWB_DAXFILE\'
  if [ -f $1 ];
  then
    if [ ! -f $CWB_DAXFILE ];
    then
      echo "output dax file : " \'$CWB_DAXFILE\' "not exist !!!"
      exit 1
    fi
  fi
fi

if [[ -z $2 ]]; then
  if [[ -z $CWB_PEGASUS_SITE ]]; then
    echo ""
    echo "Error : missing enviroment variable CWB_PEGASUS_SITE"
    echo ""
    exit 1
  else
    SITE=$CWB_PEGASUS_SITE
  fi
else 
  SITE=$2
fi

echo "DAXFILE  = " $CWB_DAXFILE 
echo "SITE     = " $SITE 

set -e

# plan and submit the  workflow
pegasus-plan \
    --conf $HOME_WAT/tools/pegasus.conf \
    --sites $SITE \
    --staging-site $SITE \
    --dir workflows \
    --output-site local \
    --dax $CWB_DAXFILE \
    -v --submit

