#!/bin/tcsh -f

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_setpars
  exit
endif

if ($1 == 'q') then
  echo ""
  echo "Current CWB_UPARAMETERS_FILE :"
  echo -n "-> "
  echo ${CWB_UPARAMETERS_FILE}
  echo ""
  exit
endif

if ($1 == '.') then
  setenv CWB_UPARAMETERS_FILE  config/user_parameters.C
else
  setenv CWB_UPARAMETERS_FILE  $1
endif

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

echo ""
echo "Current CWB_UPARAMETERS_FILE :"
echo -n "-> "
echo ${CWB_UPARAMETERS_FILE}
echo ""
