#!/bin/bash

if [[ -z $1 ]]; then
  echo ""
  echo "cwb_lsf_submit.sh :  missing dag file"
  echo ""
  exit 1
else
  CWB_DAGFILE=$1                              # DAG FILE
  echo "input  dag file : " \'$CWB_DAGFILE\'
  CWB_LSFFILE="${CWB_DAGFILE/.dag/.lsf}"      # LSF FILE
  echo "output lsf file : " \'$CWB_LSFFILE\'
  if [ -f $1 ];
  then
    if [ ! -f $CWB_LSFFILE ];
    then
      echo "output lsf file : " \'$CWB_LSFFILE\' "not exist !!!"
      exit 1
    fi
  fi
fi

echo "LSFFILE  = " $CWB_LSFFILE 

set -e

CWB_CONDOR_DIR=${PWD##*/} 
cd ..

# exit if there are jobs belonging to this group which are running
basename $PWD | awk -v user="$USER" 'BEGIN { OFS = ""; ORS = "" } ; {print "/"} ; {print user} ; {print "/"} ; {print $1} ' | xargs bjobs -g |& awk '{if($0!="No unfinished job found") {print "\nWarning : there are jobs belonging to this group which are running !!!\nBefore to submit you need to kill the running jobs : cwb_lsf kill\n"} if($0!="No unfinished job found") exit 1}' 

# submit the LSF job file
source $CWB_CONDOR_DIR/$CWB_LSFFILE

