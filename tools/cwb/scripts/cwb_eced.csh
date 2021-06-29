#!/bin/tcsh -f

onintr irq_ctrlc

if ($# == 0) then
  $CWB_SCRIPTS/cwb_help.csh cwb_eced
  exit
endif

unsetenv CWB_JOBID
unsetenv CWB_GPS_EVENT
unsetenv CWB_INET_OPTIONS
unsetenv CWB_JOB_LAG
unsetenv CWB_MDC_FACTOR
unsetenv CWB_BATCH
unsetenv CWB_CED_DIR

unsetenv CWB_ECED_OPTS
unsetenv CWB_ECED_NDET
unsetenv CWB_ECED_DET1
unsetenv CWB_ECED_DET2
unsetenv CWB_ECED_DET3
unsetenv CWB_ECED_DET4

set NDET = $#
@ NDET -= 1
if ($NDET > 0) then
  setenv CWB_ECED_NDET $NDET
endif

if ($# > 0) then
  setenv CWB_ECED_OPTS "$1"
endif

if ($NDET > 0) then
  setenv CWB_ECED_DET1 "$2"
endif
if ($NDET > 1) then
  setenv CWB_ECED_DET2 "$3"
echo $CWB_ECED_DET2
endif
if ($NDET > 2) then
  setenv CWB_ECED_DET3 "$4"
endif
if ($NDET > 3) then
  setenv CWB_ECED_DET4 "$5"
endif

setenv CWB_JOBID 1
setenv CWB_INET_OPTIONS "ced"
if ( "$1" =~ *--ced*false* ) then
  unsetenv CWB_INET_OPTIONS
endif

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_xnet.C\(\"\",0,\"\",false,true\)

unsetenv CWB_ECED_OPTS
unsetenv CWB_ECED_NDET
unsetenv CWB_ECED_DET1
unsetenv CWB_ECED_DET2
unsetenv CWB_ECED_DET3
unsetenv CWB_ECED_DET4

unsetenv CWB_JOBID
unsetenv CWB_GPS_EVENT
unsetenv CWB_INET_OPTIONS

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

