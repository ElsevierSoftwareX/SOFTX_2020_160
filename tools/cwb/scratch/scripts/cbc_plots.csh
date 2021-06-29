#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  echo ""
  echo "cb_plots mlabel bkg_tree(optional) "
  echo ""
  echo "            this command produces cbc plots + html page"
  echo "            when given the bkg tree it produces also FAD plots"
  echo ""
  echo "            mlabel : merge label"
  echo " 			bkg_tree: entire path to the the tree root file   "
   echo "            Ex :   cbc_plots M1 ../TEST1/merge/wave_TEST1.M2.root "
  echo "  "
  exit

else
  setenv CWB_MERGE_LABEL "$1"
  setenv CWB_LAG_NUMBER -1
  setenv CWB_SLAG_NUMBER -1
  if ( $2 != '' ) then
	setenv CWB_BKG_TREE "$2"
  endif
endif

  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/cbc/cbc_plots_sim4.C
 # root -n -l ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/cbc/cbc_plots_sim4.C
  if ( $? != 0) then
    echo ""
    echo "cbc_plots.C error : process terminated"
    echo ""
	exit
  endif

  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_index.C
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_index.C error : process terminated"
    echo ""
    exit
  endif

  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_header.C
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_header.C error : process terminated"
    echo ""
    exit
  endif
  
  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/cbc/cwb_mkhtml_cbc.C
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_ebbh.C error : process terminated"
    echo ""
    exit
  endif

  unsetenv CWB_MERGE_LABEL
  unsetenv CWB_BKG_TREE

exit 0
irq_ctrlc:
  ps r | grep root | awk '{print $1}' | xargs kill -9
  exit 1

