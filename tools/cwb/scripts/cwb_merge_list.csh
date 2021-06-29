#!/bin/tcsh -f

onintr irq_ctrlc

if ($1 == '') then
  echo ""
  echo "cwb_merge_list label(MX)"
  echo ""
  exit
endif

setenv CWB_MERGE_LABEL $1

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_merge_list.C

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

