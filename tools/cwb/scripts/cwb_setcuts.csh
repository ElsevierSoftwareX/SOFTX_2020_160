#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $1 '$2' $3 $4 $5 $6"

if ((($1 == '') || ("$2" == '')) && ($1 != 'list')) then
  $CWB_SCRIPTS/cwb_help.csh cwb_setcuts
  exit
endif

if ($1 == 'list') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_merge_dir.C
  exit
endif

if ( "$2" != '' ) then
  setenv CWB_SETCUTS_OPTIONS "$2"
else
  unsetenv CWB_SETCUTS_OPTIONS
endif

setenv CWB_MERGE_LABEL $1
setenv CWB_WCUTS_TREE "$2"
setenv CWB_CUTS_LABEL  $3
setenv CWB_JCUTS_TREE "$4"
setenv CWB_MCUTS_TREE "$5"
setenv CWB_LCUTS_TREE "$6"
setenv CWB_LAG_NUMBER -1
setenv CWB_SLAG_NUMBER -1

root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_setcuts.C
if ( $? != 0) exit 1

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

unsetenv CWB_MERGE_LABEL
unsetenv CWB_WCUTS_TREE
unsetenv CWB_CUTS_LABEL
unsetenv CWB_JCUTS_TREE
unsetenv CWB_MCUTS_TREE
unsetenv CWB_LCUTS_TREE
unsetenv CWB_SETCUTS_OPTIONS
unsetenv CWB_SLAG_NUMBER
unsetenv CWB_LAG_NUMBER

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

