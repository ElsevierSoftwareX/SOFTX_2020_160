#!/bin/tcsh -f

onintr irq_ctrlc

 setenv CWB_XPARAMETERS_FILE ${CWB_UPARAMETERS_FILE}

 setenv CWB_XPPARAMETERS_FILE ${CWB_UPPARAMETERS_FILE}

set CWB_XPARM  = "setenv CWB_UPARAMETERS_FILE  ${CWB_XPARAMETERS_FILE}"
set CWB_XPPARM = "setenv CWB_UPPARAMETERS_FILE ${CWB_XPPARAMETERS_FILE}"

setenv CWB_MERGE_LABEL M1
 setenv OUTPUT_DIR "$1"

  (${CWB_XPARM};${CWB_XPPARM};root -n -l ${CWB_PARMS_FILES} ${CWB_MACROS}/../postproduction/burst/my_condor_benchmark.C)


exit 0
irq_ctrlc:
  ps r | grep root | awk '{print $1}' | xargs kill -9
  exit 1

