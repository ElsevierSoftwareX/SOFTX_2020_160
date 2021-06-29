#!/bin/bash

source $7

export CWB_JOBID=$1
export CWB_UFILE=$2
export CWB_STAGE=$3
export CWB_WDIR=$4
export CWB_TGZFILE=$5
export CWB_OUTFILE=$6

pwd 

export _USE_PEGASUS=1

/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh $CWB_WDIR
cd $CWB_WDIR
tar -xzvf ../$CWB_TGZFILE > /dev/null
mv ../*.root output/
tar cvf /dev/null . > CWB_TGZFILE

echo ${ROOTSYS}/bin/root -b -q -l ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_DIR_NAME}/cwb_xnet.C\(\"${CWB_UFILE}\",${CWB_STAGE},\"${CWB_UPARAMETERS_FILE}\",true,false\)

${ROOTSYS}/bin/root -b -q -l ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_DIR_NAME}/cwb_xnet.C\(\"${CWB_UFILE}\",${CWB_STAGE},\"${CWB_UPARAMETERS_FILE}\",true,false\)

CWB_OUTTGZ=$CWB_OUTFILE".job"$CWB_JOBID".tgz"
tar czvf $CWB_OUTTGZ -X CWB_TGZFILE log/ output/

ls -laR

mv $CWB_OUTTGZ ../
cd ..

