#!/bin/tcsh -f

onintr irq_ctrlc

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_mkfad
  exit
endif

setenv CWB_MERGE_LABEL $1
setenv CWB_LAG_NUMBER -1
setenv CWB_SLAG_NUMBER -1

if ($2 == '') then
  setenv CWB_MKFAD_NZBINS -1
else
  if ( $2 =~ *[^0-9]* ) then
    setenv CWB_FAD_CONFIG $2
    setenv CWB_MKFAD_NZBINS -1
  else
    setenv CWB_MKFAD_NZBINS $2
    setenv CWB_FAD_CONFIG ""
  endif
endif

if ($3 == '') then
  setenv CWB_MKFAD_DIR "fad"
else
  setenv CWB_MKFAD_DIR $3
endif

if ($4 == '') then
  setenv CWB_MKFAD_LIVETIME 0
else
  setenv CWB_MKFAD_LIVETIME $4
endif

root -n -l -b ${CWB_PPARMS_FILES} ${CWB_FAD_CONFIG} ${CWB_MACROS}/cwb_mkfad.C\(\"$CWB_MKFAD_DIR\",$CWB_MKFAD_NZBINS,$CWB_MKFAD_LIVETIME\)

unsetenv CWB_MERGE_LABEL
unsetenv CWB_LAG_NUMBER 
unsetenv CWB_SLAG_NUMBER 
unsetenv CWB_MKFAD_NZBINS
unsetenv CWB_MKFAD_DIR
unsetenv CWB_MKFAD_LIVETIME
unsetenv CWB_FAD_CONFIG

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

