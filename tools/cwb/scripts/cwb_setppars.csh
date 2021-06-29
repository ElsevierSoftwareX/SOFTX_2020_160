#!/bin/tcsh -f

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_setppars
  exit
endif

if ($1 == 'q') then
  echo ""
  echo "Current CWB_UPPARAMETERS_FILE :"
  echo -n "-> "
  echo ${CWB_UPPARAMETERS_FILE}
  echo ""
  exit
endif

if ($1 == '.') then
  setenv CWB_UPPARAMETERS_FILE  config/user_pparameters.C
else
  setenv CWB_UPPARAMETERS_FILE  $1
endif

echo ""
echo "Current CWB_UPPARAMETERS_FILE :"
echo -n "-> "
echo ${CWB_UPPARAMETERS_FILE}
echo ""
