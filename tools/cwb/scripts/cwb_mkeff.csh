#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  echo "cwb_mkeff post-prod-label (ex: postprod.M1)"
  exit
endif

root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkeff.C

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

