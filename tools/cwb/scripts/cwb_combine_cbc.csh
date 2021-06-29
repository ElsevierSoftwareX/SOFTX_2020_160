#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $1 '$2'"

if ((($1 == '') || ("$2" == '')) && ($1 != 'list')) then
#  $CWB_SCRIPTS/cwb_help.csh cwb_combine_cbc
  echo ""
  echo "Available options:"
  echo ""
  echo " --fthr    -> frequency threshold used to combine IMBHB & BBH searches"
  echo " --search  -> search type (IMBHB or BBH) to be combined "
  echo " --wdir    -> the working dir of the search to be combined"
  echo " --mlabel  -> the merge label of the search to be combined"
  echo " --ulabel  -> the user label to be used to tag the combined output wave/mdc files"
  echo " --lag     -> lag used to select background data to be combined"
  echo " --slag    -> slag used to select background data to be combined"
  echo " --ifarthr -> only for background. Used to select events with ifar>ifarthr, ifarthr is in years"
  echo ""
  echo "Ex: cwb_combine_cbc M1.C_U.S_bin1_cut '--search IMBHB --wdir /home/IMBHB/LH/SIM/O3_K02_C00_LH_IMBHB_SIM_NR-MIX_run1 --mlabel M1.C_U.S_bin1_cut --fthr 80'"
  echo ""
  exit
endif

if ($1 == 'list') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_merge_dir.C
  exit
endif

if ( "$2" != '' ) then
  setenv CWB_COMBINE_OPTIONS "$2"
else
  unsetenv CWB_COMBINE_OPTIONS
endif

setenv CWB_MERGE_LABEL $1

root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_combine_cbc.C
if ( $? != 0) exit 1

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

unsetenv CWB_MERGE_LABEL

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

