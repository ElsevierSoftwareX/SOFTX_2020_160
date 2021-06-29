#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  $CWB_SCRIPTS/cwb_help.csh cwb_draw_sensitivity
  exit
endif

setenv CWB_SENSITIVITY_FILE_NAME $1
setenv CWB_SENSITIVITY_SAVE_PLOT $2
setenv CWB_SENSITIVITY_RANGE_FIX $3
root -n -l ${CWB_MACROS}/cwb_draw_sensitivity.C 
unsetenv CWB_SENSITIVITY_FILE_NAME 
unsetenv CWB_SENSITIVITY_SAVE_PLOT 
unsetenv CWB_SENSITIVITY_RANGE_FIX 

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

