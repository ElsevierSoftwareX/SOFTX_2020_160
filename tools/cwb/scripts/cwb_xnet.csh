#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  echo "cwb_inet jobid gps_evt[0=(full job)/sec] options[false/true] mdc_factor lag[0=(full lags)/Xlag] ced_dir(def:data)"
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
  setenv CWB_GPS_EVENT $2
else
  unsetenv CWB_GPS_EVENT
endif

if ( "$3" != '' ) then
  setenv CWB_INET_OPTIONS "$3"
else
  unsetenv CWB_INET_OPTIONS
endif

if ( $4 != '' ) then
  setenv CWB_MDC_FACTOR $4
else
  unsetenv CWB_MDC_FACTOR 
endif

if ( $5 != '' ) then
  setenv CWB_JOB_LAG $5
else
  unsetenv CWB_JOB_LAG 
endif

if ( $6 != '' ) then
  setenv CWB_CED_DIR $6
else
  unsetenv CWB_CED_DIR 
endif

root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_MACROS}/cwb_xnet.C\(\"${CWB_UPARAMETERS_FILE}\"\)

unsetenv CWB_JOBID
unsetenv CWB_GPS_EVENT
unsetenv CWB_INET_OPTIONS 
unsetenv CWB_JOB_LAG 
unsetenv CWB_MDC_FACTOR 
unsetenv CWB_CED_DIR 
unsetenv CWB_DUMP_INFOS_AND_EXIT

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

