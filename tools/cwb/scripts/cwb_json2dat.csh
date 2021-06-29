#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_json2dat
  exit
endif

if ( ! -d $1 ) then 
  echo ""
  echo "file $1 not exist"
  echo ""
  exit
endif

python $CWB_SCRIPTS/cwb_json2dat.py $1
if ( $? != 0) exit 1

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

