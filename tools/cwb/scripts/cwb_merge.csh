#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_merge
  exit
endif

if ($1 == 'list') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_merge_dir.C
  exit
endif

setenv CWB_MERGE_LABEL $1
setenv CWB_MERGE_OPTS "$2"

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_merge.C
if ( $? != 0) exit 1
# makes ROOT5 macro configuration files ROOT6 compliant 
root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_MACROS}/FixUserParametersROOT6.C\(\"$CWB_MPARAMETERS_FILE\",\"mp\"\)

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

unsetenv CWB_MERGE_LABEL
unsetenv CWB_MERGE_OPTS

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

