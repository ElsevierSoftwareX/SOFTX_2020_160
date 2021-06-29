#!/bin/bash

cd ..
${ROOTSYS}/bin/root -b -q -l ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_xnet.C\(\"${CWB_UFILE}\",${CWB_STAGE},\"${CWB_UPARAMETERS_FILE}\",true,false\)


#exit 0
