#!/bin/tcsh -f

onintr irq_ctrlc

if ($1 == '') then
  echo ""
  echo "cwb_fix_slag_missed label/list (ex: M1)"
  echo ""
  exit
endif

if ($1 == 'list') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_merge_dir.C
  exit
endif

setenv CWB_MERGE_LABEL $1

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_fix_live_slag_missed.C

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

