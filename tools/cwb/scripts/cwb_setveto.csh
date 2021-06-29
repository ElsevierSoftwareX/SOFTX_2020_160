#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_setveto
  exit
endif

setenv CWB_LAG_NUMBER -1
setenv CWB_SLAG_NUMBER -1

if ($1 == 'list') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_merge_dir.C
  exit
endif

setenv CWB_MERGE_LABEL $1

if ($2 != '') then
  setenv CWB_XPARAMETERS_FILE $2
else
  setenv CWB_XPARAMETERS_FILE ${CWB_UPARAMETERS_FILE}
endif

if ($3 != '') then
  setenv CWB_XPPARAMETERS_FILE $3
else
  setenv CWB_XPPARAMETERS_FILE ${CWB_UPPARAMETERS_FILE}
endif

set CWB_XPARM  = "setenv CWB_UPARAMETERS_FILE  ${CWB_XPARAMETERS_FILE}"
set CWB_XPPARM = "setenv CWB_UPPARAMETERS_FILE ${CWB_XPPARAMETERS_FILE}"

(${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_setveto.C)
if ( $? != 0) exit 1

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

unsetenv CWB_MERGE_LABEL
unsetenv CWB_LAG_NUMBER
unsetenv CWB_SLAG_NUMBER
unsetenv CWB_XPARAMETERS_FILE
unsetenv CWB_XPPARAMETERS_FILE

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

