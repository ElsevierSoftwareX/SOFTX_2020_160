#!/bin/tcsh -f

if (($1 != '1G') && ($1 != '1g') && ($1 != '1') && ($1 != '2G') && ($1 != '2g') && ($1 != '2')) then
  $CWB_SCRIPTS/cwb_help.csh cwb_setpipe
  echo ""
  echo -n "Current CWB_ANALYSIS="
  echo ${CWB_ANALYSIS}
  echo ""
  exit
endif

if (($1 == '1G') || ($1 == '1g') || ($1 == '1')) then
  setenv CWB_ANALYSIS        "1G"       # 1G  analysis
  source $HOME_WAT/tools/config.csh
  echo ""
  echo -n "CWB_ANALYSIS="
  echo ${CWB_ANALYSIS}
  echo ""
endif

if (($1 == '2G') || ($1 == '2g') || ($1 == '2')) then
  setenv CWB_ANALYSIS        "2G"       # 2G  analysis
  source $HOME_WAT/tools/config.csh
  echo ""
  echo -n "CWB_ANALYSIS="
  echo ${CWB_ANALYSIS}
  echo ""
endif

