#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  echo "cwb_inet jobid/fName data_dir(def:data)"
  exit
endif

unsetenv CWB_JOBID
unsetenv CWB_GPS_EVENT
unsetenv CWB_INET_OPTIONS
unsetenv CWB_JOB_LAG
unsetenv CWB_MDC_FACTOR
unsetenv CWB_BATCH
unsetenv CWB_CED_DIR

setenv CWB_JOBID $1

if ( $2 != '' ) then
  setenv CWB_CED_DIR $2
else
  unsetenv CWB_CED_DIR 
endif

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_xnet.C\(\"$CWB_JOBID\",CWB_STAGE_LIKELIHOOD,\"${CWB_UPARAMETERS_FILE}\"\)

unsetenv CWB_JOBID
unsetenv CWB_CED_DIR 

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

