#!/bin/tcsh -f

onintr irq_ctrlc

if (($1 == '') || ($2 == '')) then
  echo "cwb_frdisplay jobid ifo(H1,L1,V1,G1)"
  exit
endif
setenv CWB_JOBID $1
setenv CWB_IFO $2
root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_frdisplay.C 

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

