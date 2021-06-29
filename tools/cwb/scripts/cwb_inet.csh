#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  $CWB_SCRIPTS/cwb_help.csh cwb_inet
  exit
endif

unsetenv CWB_JOBID
unsetenv CWB_GPS_EVENT
unsetenv CWB_INET_OPTIONS 
unsetenv CWB_JOB_LAG 
unsetenv CWB_MDC_FACTOR 
unsetenv CWB_UFILE_NAME 
unsetenv CWB_BATCH 
unsetenv CWB_CED_DIR 

setenv CWB_JOBID $1
setenv CWB_UFILE_NAME ${CWB_UPARAMETERS_FILE}

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
  ## check if $6 is a directory than it is used as ced dir
  ## otherwise it used as the input user_parameters.C file
  if (! -d $6) then
    if ( ! -e $6 ) then
      echo ""
      echo -n "configuration file '" $6
      echo " ' not exist!"
      echo ""
      exit
    else
      setenv CWB_UFILE_NAME $6
    endif 
  else
    setenv CWB_CED_DIR $6
  endif
else
  unsetenv CWB_CED_DIR 
endif

if (! $?_USE_ROOT6 ) then
root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_xnet.C\(\"${CWB_UFILE_NAME}\"\)
else
root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_MACROS}/cwb_xnet.C\(\"${CWB_UFILE_NAME}\"\)
endif

unsetenv CWB_JOBID
unsetenv CWB_GPS_EVENT
unsetenv CWB_INET_OPTIONS 
unsetenv CWB_JOB_LAG 
unsetenv CWB_MDC_FACTOR 
unsetenv CWB_UFILE_NAME 
unsetenv CWB_CED_DIR 

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

