#!/bin/tcsh -f

onintr irq_ctrlc

if (($1 == '') || (( $1 != 'create' ) && ( $1 != 'log' ))) then
  echo "wmdc_condor [create/log]"
  exit
endif

if ( $1 == 'create' ) then

  root -n -l -b ${CWB_ROOTLOGON_FILE} ${WMDC_CONFIG} ${HOME_WAVEMDC}/wmdc_condor_create.C

endif

if ( $1 == 'log' ) then

  root -n -l -b ${CWB_ROOTLOGON_FILE} ${WMDC_CONFIG} ${HOME_WAVEMDC}/wmdc_condor_log.C

endif

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

