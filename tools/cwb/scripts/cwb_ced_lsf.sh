#!/bin/bash

export CWB_JOBID=$1
export CWB_STAGE=$2
export CWB_NODEDIR=$3
export CWB_DATA_LABEL=$4
export CWB_DATA_DIR=$5
export CWB_UFILE=$6

# make lsf label
if [[ $CWB_STAGE == "CWB_STAGE_FULL" ]]; then
  export CWB_LSF_LABEL=$CWB_JOBID\_$CWB_DATA_LABEL
else
  export CWB_LSF_LABEL=$CWB_JOBID\_$CWB_DATA_LABEL\_$CWB_STAGE
fi

# cleanup nodedir
rm -rf $CWB_NODEDIR/$CWB_LSF_LABEL
cd $CWB_NODEDIR

# create working dir
mkdir $CWB_LSF_LABEL
cd $CWB_LSF_LABEL

# copy working directories from tgz
tar -xzf ../$CWB_LSF_LABEL.tgz

# remove tgz
rm ../$CWB_LSF_LABEL.tgz

# create data dir
mkdir $CWB_DATA_DIR

# copy ufile 
mkdir -p `dirname $CWB_UFILE`
cp ../$CWB_LSF_LABEL.ufile $CWB_UFILE

# remove original ufile
rm ../$CWB_LSF_LABEL.ufile

# define pipeline command
CWB_NET="${ROOTSYS}/bin/root -b -q -l ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_xnet.C(\"${CWB_UPARAMETERS_FILE}\",0,"false",false)"

# execute analysis
echo $CWB_NET
$CWB_NET

# produce tgz of the analysis output 
cd $CWB_DATA_DIR
tar -czf $CWB_NODEDIR/$CWB_LSF_LABEL.tgz *
cd ../..

# remove working directory
rm -rf $CWB_NODEDIR/$CWB_LSF_LABEL
