#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  echo "please provide config file" 
  echo "Ex: cwb_cmpeff config/cwb_cmpeff_config.C"
  exit
endif

root -n -l -b ${CWB_PPARMS_FILES} $1 ${CWB_MACROS}/cwb_cmpeff.C

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

