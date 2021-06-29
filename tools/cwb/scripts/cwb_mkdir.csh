#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_mkdir
  exit
endif

if (($2 != '') && ($2 != 'batch') && ($2 != 'interactive')) then
  echo "cwb_mkdir wrk_dir option [interactive/batch]"
  exit
endif

unset CWB_MKDIR_WRKDIR
if ($1 == '.') then
  setenv CWB_MKDIR_WRKDIR ""
else
  setenv CWB_MKDIR_WRKDIR $1
endif
unset CWB_MKDIR_OPTION
if ($2 == '') then
  setenv CWB_MKDIR_OPTION "batch"
else
  setenv CWB_MKDIR_OPTION $2
endif

root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_MACROS}/cwb_mkdir.C
if ( $? != 0) exit 1

# create cWB_analysis.log file
cd $CWB_MKDIR_WRKDIR
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

echo ''
echo 'The new working dir is :  '$CWB_MKDIR_WRKDIR
echo ''

cd ..

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

