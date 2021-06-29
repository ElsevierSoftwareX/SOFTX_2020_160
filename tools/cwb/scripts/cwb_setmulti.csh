#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  echo ""
  echo "cwb_setmulti label/list (ex: M1)"
  echo ""
  exit
endif

if ($1 == 'list') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_merge_dir.C
  exit
endif

setenv CWB_MERGE_LABEL $1

root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_setmulti.C
if ( $? != 0) exit 1

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

