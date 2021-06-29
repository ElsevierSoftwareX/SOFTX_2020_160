#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' || $2 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_mplugin
  exit
endif

unsetenv CWB_MPLUGIN_1
unsetenv CWB_MPLUGIN_2
unsetenv CWB_MPLUGIN_3
unsetenv CWB_MPLUGIN_4
unsetenv CWB_MPLUGIN_5
unsetenv CWB_MPLUGIN_6

if ( $1 != '' ) then
  setenv CWB_MPLUGIN_1 $1
endif
if ( $2 != '' ) then
  setenv CWB_MPLUGIN_2 $2
endif
if ( $3 != '' ) then
  setenv CWB_MPLUGIN_3 $3
endif
if ( $4 != '' ) then
  setenv CWB_MPLUGIN_4 $4
endif
if ( $5 != '' ) then
  setenv CWB_MPLUGIN_5 $5
endif
if ( $6 != '' ) then
  setenv CWB_MPLUGIN_6 $6
endif
if ( $7 != '' ) then
  echo ""
  echo "cwb_mplugin error : max number of merged plugins are 5"
  echo ""
  exit
endif

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_mplugin.C

unsetenv CWB_MPLUGIN_1
unsetenv CWB_MPLUGIN_2
unsetenv CWB_MPLUGIN_3
unsetenv CWB_MPLUGIN_4
unsetenv CWB_MPLUGIN_5
unsetenv CWB_MPLUGIN_6

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

