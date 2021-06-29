#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  echo "cwb_jnet jobFile jobConf[optional]"
  exit
endif

if ( $1 != '' ) then
  setenv CWB_JOB_FILE $1
else
  setenv CWB_JOB_FILE "" 
endif

if ( $2 != '' ) then
  setenv CWB_JOB_CONF $2
else
  setenv CWB_JOB_CONF ""
endif

root -n -l -b  ${CWB_ROOTLOGON_FILE} ${CWB_MACROS}/cwb_jnet.C\(\"${CWB_JOB_FILE}\"\,\"${CWB_JOB_CONF}\"\)

unsetenv CWB_JOB_FILE
unsetenv CWB_JOB_CONF

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

