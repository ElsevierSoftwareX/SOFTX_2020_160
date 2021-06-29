#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_petools
  exit
endif

if ($1 == '') then
  setenv CWB_IFILE ""
else
  setenv CWB_IFILE $1
endif

if ($2 == '') then
  setenv CWB_OPTION ""
else
  setenv CWB_OPTION $2
endif

if ($3 == '') then
  setenv CWB_OPATH ""
else
  setenv CWB_OPATH $3
endif

if ($4 == '') then
  setenv CWB_PE_GIT "/home/waveburst/git/pe"
else
  setenv CWB_PE_GIT $4
endif

root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_MACROS}/cwb_petools.C\(\"$CWB_IFILE\",\"$CWB_OPTION\",\"$CWB_OPATH\",\"$CWB_PE_GIT\"\)
if ( $? != 0) exit 1

# create cWB_analysis.log file
#make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

unsetenv CWB_IFILE
unsetenv CWB_OPTION
unsetenv CWB_OPATH
unsetenv CWB_PE_GIT

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

