#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ((($1 == '') || ("$2" == '')) && ($1 != 'list')) then
  $CWB_SCRIPTS/cwb_help.csh cwb_setifar
  exit
endif

if ($1 == 'list') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_merge_dir.C
  exit
endif

if ( "$2" != '' ) then
  setenv CWB_SETIFAR_OPTIONS "$2"
else
  unsetenv CWB_SETIFAR_OPTIONS
endif

setenv CWB_MERGE_LABEL     $1
setenv CWB_SETIFAR_TSEL   "$2"
setenv CWB_SETIFAR_FILE    $3
setenv CWB_SETIFAR_LABEL   $4

if ( "$5" != '' ) then
  setenv CWB_SETIFAR_MODE "$5"
else
  unsetenv CWB_SETIFAR_MODE
endif

root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_setifar.C
if ( $? != 0) exit 1

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

unsetenv CWB_MERGE_LABEL
unsetenv CWB_SETIFAR_TSEL
unsetenv CWB_SETIFAR_FILE
unsetenv CWB_SETIFAR_LABEL
unsetenv CWB_SETIFAR_MODE
unsetenv CWB_SETIFAR_OPTIONS

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

