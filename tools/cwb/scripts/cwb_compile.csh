#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_compile
  exit
endif

setenv CWB_COMPILE_MACRO_NAME $1

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_compile.C

unsetenv CWB_COMPILE_MACRO_NAME

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

