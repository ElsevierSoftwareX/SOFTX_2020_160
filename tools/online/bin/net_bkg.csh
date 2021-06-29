#!/bin/tcsh -f

setenv CWB_UPARAMETERS_FILE ${CWB_UFILE}
setenv Slag_datashift ${SLAG_SHIFT}
#${HOME_CWB}/scripts/cwb_net.sh
#${HOME_CWB}/scripts/cwb_inet.csh ${CWB_JOBID}
cd ..
${ROOTSYS}/bin/root -b -q -l -n ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_xnet.C\(\"${CWB_UFILE}\",${CWB_STAGE},\"${CWB_UPARAMETERS_FILE}\",true,false\)

