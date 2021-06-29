#!/bin/tcsh -f

onintr irq_ctrlc

if (($1 == '') || ($2 == '')) then
  $CWB_SCRIPTS/cwb_help.csh cwb_xtalk
  exit
endif


setenv CWB_XTALK_LOW_LEVEL  $1
setenv CWB_XTALK_HIGH_LEVEL $2

if ($3 != '') then
  setenv CWB_XTALK_INU $3
else
unsetenv CWB_XTALK_INU
endif

if ($4 != '') then
  setenv CWB_XTALK_PRECISION $4
else
unsetenv CWB_XTALK_PRECISION
endif

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_xtalk.C

unsetenv CWB_XTALK_LOW_LEVEL
unsetenv CWB_XTALK_HIGH_LEVEL
unsetenv CWB_XTALK_INU
unsetenv CWB_XTALK_PRECISION

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

