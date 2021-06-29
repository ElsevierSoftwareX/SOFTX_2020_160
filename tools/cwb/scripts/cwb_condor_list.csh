#!/bin/tcsh -f

onintr irq_ctrlc

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_list.C

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

