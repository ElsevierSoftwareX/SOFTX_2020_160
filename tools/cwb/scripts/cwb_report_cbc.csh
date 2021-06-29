#!/bin/tcsh -f

onintr irq_ctrlc

if ( $1 == '' ) then
  echo ""
  echo "cb_plots_ifar mlabel "
  echo ""
  echo "            this command produces cbc plots + html page for sim4 simulations"
  echo "            when given a tree with ifar branch it produces also FAD plots"
  echo ""
  echo "            mlabel : merge label"
#  echo " 			bkg_tree: entire path to the the tree root file   "
   echo "            Ex :   cbc_plots_ifar M1 "
  echo "  "
  exit

else
  setenv CWB_MERGE_LABEL "$1"
  setenv CWB_LAG_NUMBER -1
  setenv CWB_SLAG_NUMBER -1
  setenv CWB_REPORT_CBC 1
  if ( $2 != '' ) then
	setenv CBC_SEARCH_TYPE "$2"
  else      
    setenv CBC_SEARCH_TYPE "BBH"    
  endif
endif


  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_report_cbc.C
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
  
  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_cbc.C
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_cbc.C error : process terminated"
    echo ""
    exit
  endif

  unsetenv CWB_MERGE_LABEL
  unsetenv CWB_LAG_NUMBER 
  unsetenv CWB_SLAG_NUMBER 
  unsetenv CWB_REPORT_CBC 
  
exit 0
irq_ctrlc:
  ps r | grep root | awk '{print $1}' | xargs kill -9
  exit 1

