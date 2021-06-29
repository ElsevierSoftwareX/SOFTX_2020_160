#!/bin/bash

${ROOTSYS}/bin/root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ./macro/Init.C ./macro/CWB_Plugin_OnlineFrame_Config.C macro/mlog.C
