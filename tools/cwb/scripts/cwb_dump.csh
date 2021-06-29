#!/bin/tcsh -f

onintr irq_ctrlc

if (($1 == '') || (($1 != 'all') && ($1 != 'dq') && ($1 != 'inj') && ($1 != 'job') && ($1 != 'sjob') && ($1 != 'lag') && ($1 != 'slag') && ($1 != 'history') && ($1 != 'events') && ($1 != 'config'))) then
  $CWB_SCRIPTS/cwb_help.csh cwb_dump
  exit
endif

if ( $1 == 'all' ) then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_dq.C
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_inj.C
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_lag.C
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_slag.C
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_job.C
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_sjob.C
endif

if ( $1 == 'dq' ) then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_dq.C
endif

if ( $1 == 'inj' ) then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_inj.C
endif

if ( $1 == 'lag' ) then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_lag.C
endif

if ( $1 == 'slag' ) then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_slag.C
endif

if ( $1 == 'job' ) then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_job.C
endif

if ( $1 == 'sjob' ) then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_sjob.C
endif

if ( $1 == 'history' ) then

  if (($2 == '') || (($3 != 'dump') && ($3 != 'view') && ($3 != ''))) then
    echo ""
    echo "cwb_dump history file_name dump/view(def: view)"
    echo ""
    exit
  endif

  setenv CWB_DUMP_HIST_FILE_NAME $2
  unsetenv CWB_DUMP_HIST_MODE
  if ($3 != '') then
    setenv CWB_DUMP_HIST_MODE $3
  endif

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_history.C

  unsetenv CWB_DUMP_HIST_MODE
  unsetenv CWB_DUMP_HIST_FILE_NAME

endif

if ( $1 == 'config' ) then

  if (($2 == '') || (($3 != 'dump') && ($3 != 'view') && ($3 != 'md5') && ($3 != ''))) then
    echo ""
    echo "cwb_dump config file_name dump/view/md5(def: view)"
    echo ""
    echo "   - get the cwb parameters.C used in the last stage"
    echo "     - dump : save config to file report/dump/config_*.C"
    echo "     - view : view config"
    echo "     - md5  : view the MD5 cwb parameters.C"
    echo ""
    exit
  endif

  setenv CWB_DUMP_HIST_FILE_NAME $2
  unsetenv CWB_DUMP_HIST_MODE
  if ($3 != '') then
    setenv CWB_DUMP_HIST_MODE $3
  endif

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_config.C

  unsetenv CWB_DUMP_HIST_MODE
  unsetenv CWB_DUMP_HIST_FILE_NAME

endif

if ( $1 == 'events' ) then

  if (($2 == '') || (($3 != 'dump') && ($3 != 'view') && ($3 != ''))) then
    echo ""
    echo "cwb_dump events file_name dump/view(def: view)"
    echo ""
    exit
  endif

  setenv CWB_DUMP_EVT_FILE_NAME $2
  unsetenv CWB_DUMP_EVT_MODE
  if ($3 != '') then
    setenv CWB_DUMP_EVT_MODE $3
  endif

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_events.C

  unsetenv CWB_DUMP_EVT_MODE
  unsetenv CWB_DUMP_EVT_FILE_NAME

endif

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

