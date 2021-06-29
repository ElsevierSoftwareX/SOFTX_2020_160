#!/bin/tcsh -f

setenv CWB_UPARAMETERS_FILE $1
setenv Slag_datashift $2
#${HOME_CWB}/scripts/cwb_net.sh
${HOME_CWB}/scripts/cwb_inet.csh $3 $4 true 0 $5 'ced' >& log_ced_$4_$5 &

