#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  $CWB_SCRIPTS/cwb_help.csh cwb_draw_antpat
  exit
endif

if (($1 != 0) && ($1 != 1) && ($1 != 2) && ($1 != 3)) then
  setenv CWB_ANTPAT_POLARIZATION 3
else 
  setenv CWB_ANTPAT_POLARIZATION $1
endif

if (($2 != 'hammer') && ($2 != 'rect') && ($2 != '')) then
  $CWB_SCRIPTS/cwb_help.csh cwb_draw_antpat
  exit
endif

if (($2 != 'hammer') && ($2 != 'rect')) then
  setenv CWB_ANTPAT_PROJECTION "hammer"
else 
  setenv CWB_ANTPAT_PROJECTION $2
endif

if (($3 != 0) && ($3 != 1)) then
  setenv CWB_ANTPAT_SAVE_PLOT 0
else 
  setenv CWB_ANTPAT_SAVE_PLOT $3
endif

root -n -l ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_draw_antpat.C
unsetenv CWB_ANTPAT_PROJECTION 
unsetenv CWB_ANTPAT_SAVE_PLOT
unsetenv CWB_ANTPAT_POLARIZATION 

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

