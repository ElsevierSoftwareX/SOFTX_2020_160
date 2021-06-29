#!/bin/tcsh -f

onintr irq_ctrlc

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_mkrep
  exit
endif

setenv CWB_MKREP_INDEX_FILE $1
setenv CWB_LAG_NUMBER -1
setenv CWB_SLAG_NUMBER -1
setenv CWB_MERGE_LABEL ""

root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_index.C

unsetenv CWB_MKREP_INDEX_FILE
unsetenv CWB_LAG_NUMBER
unsetenv CWB_SLAG_NUMBER
unsetenv CWB_MERGE_LABEL

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

