#!/bin/tcsh -f

onintr irq_ctrlc

if (($1 == '') || (( $2 != '' ) && ( $2 != 'clean' ))) then
  echo ""
  echo "cwb_publish directory/list option [clean]"
  echo ""
  exit
endif

if ($1 == 'list') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_report_dir.C
  exit
endif

setenv CWB_PP_DATA_DIR $1
setenv CWB_PUBLISH_OPTION $2
root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_publish.C
root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_report_dir.C

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

