#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  echo "cwb_simplot merge_version (ex: M1)"
  exit
endif

setenv CWB_MERGE_LABEL $1
setenv CWB_LAG_NUMBER -1
setenv CWB_SLAG_NUMBER -1

unsetenv CWB_MKREP_INDEX_FILE
root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_index.C
if ( $? != 0) then
  echo ""
  echo "cwb_mkhtml_index.C error : process terminated"
  echo ""
  exit
endif
root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_header.C
if ( $? != 0) then
  echo ""
  echo "cwb_mkhtml_header.C error : process terminated"
  echo ""
  exit
endif
root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_sim.C
if ( $? != 0) then
  echo ""
  echo "cwb_mkhtml_sim.C error : process terminated"
  echo ""
  exit
endif

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

