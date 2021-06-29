#!/bin/bash

env
cd ..
#${ROOTSYS}/bin/root -b -q -l ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_inet.C ${CWB_NETC_FILE}
${ROOTSYS}/bin/root -b -q -l ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_xnet.C\(\"${CWB_UPARAMETERS_FILE}\",CWB_STAGE_FULL,\"\",false\)

#exit(0);
