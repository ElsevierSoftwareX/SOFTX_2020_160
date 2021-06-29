#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' || $2 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_inet2G
  exit
endif

if (($2 == '') || (($2 != 'FULL') && ($2 != 'INIT') && ($2 != 'STRAIN') && ($2 != 'CSTRAIN') && ($2 != 'COHERENCE') && ($2 != 'SUPERCLUSTER') && ($2 != 'LIKELIHOOD'))) then
  $CWB_SCRIPTS/cwb_help.csh cwb_inet2G
  exit
endif

if ( $2 != '' ) then
  setenv CWB_STAGE_NAME "CWB_STAGE_"$2
else
  setenv CWB_STAGE_NAME "CWB_STAGE_FULL"
endif

setenv CWB_UFILE_NAME ""
setenv CWB_JOBID 0
if ( $3 != '' ) then
  if ( $3 =~ *[^0-9]* ) then 
    setenv CWB_UFILE_NAME $3
  else
    setenv CWB_JOBID $3
  endif
endif

if ( "$4" == 'help' ) then
  $CWB_SCRIPTS/cwb_help.csh cwb_inet2G
  exit
endif

if ( "$4" != '' ) then
  setenv CWB_INET_OPTIONS "$4"
else
  unsetenv CWB_INET_OPTIONS
endif

if ( $5 != '' ) then
  setenv CWB_MDC_FACTOR $5
else
  unsetenv CWB_MDC_FACTOR
endif

if ( $6 != '' ) then
  setenv CWB_JOB_LAG $6
else
  unsetenv CWB_JOB_LAG
endif

if ( $7 != '' ) then
  ## check if $7 is a directory than it is used as ced dir
  ## otherwise it used as the input user_parameters.C file
  if (! -d $7) then
    if ( ! -e $7 ) then
      echo ""
      echo -n "configuration file '" $7
      echo " ' not exist!"
      echo ""
      exit
    else
      setenv CWB_UFILE_NAME $7
    endif
  else
    setenv CWB_CED_DIR $7
  endif
else
  unsetenv CWB_CED_DIR
endif

if (! $?_USE_ROOT6 ) then
root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_MACROS}/cwb_inet2G.C\($CWB_JOBID,\"$1\",$CWB_STAGE_NAME,\"$CWB_UFILE_NAME\"\)
else
root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_MACROS}/cwb_inet2G.C\($CWB_JOBID,\"$1\",$CWB_STAGE_NAME,\"$CWB_UFILE_NAME\"\)
endif

unsetenv CWB_JOBID
unsetenv CWB_STAGE_NAME
unsetenv CWB_UFILE_NAME
unsetenv CWB_INET_OPTIONS
unsetenv CWB_JOB_LAG
unsetenv CWB_MDC_FACTOR
unsetenv CWB_CED_DIR

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

