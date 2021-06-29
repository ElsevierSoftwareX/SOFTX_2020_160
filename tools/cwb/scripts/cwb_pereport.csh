#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_pereport
  exit
endif

make -f $CWB_SCRIPTS/Makefile.cwb_pereport $1 $2 $3 $4 $5 $6 $8 $9 $10 $11
if ( $? != 0) exit 1

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

