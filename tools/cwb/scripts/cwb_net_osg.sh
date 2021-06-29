#!/bin/bash

export WORKING_DIR="$(ls *.tgz)"
export WORKING_DIR="$(basename $WORKING_DIR .tgz)"
mkdir $WORKING_DIR
export HOME="${PWD}/${WORKING_DIR}"
mv *.tgz $HOME/
mv .rootrc $HOME/
cd $HOME/

tar -zxf *.tgz

${ROOTSYS}/bin/root -b -q -l ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_xnet.C\(\"${CWB_UFILE}\",${CWB_STAGE},\"${CWB_UPARAMETERS_FILE}\",false,false\)


